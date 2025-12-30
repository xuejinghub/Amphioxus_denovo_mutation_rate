#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys, csv, re, argparse
from collections import defaultdict
from typing import Optional, List, Tuple, Dict

BUBBLE_PAT   = re.compile(r'^(primary|secondary)_bubble(\d+).*?([MP])$')
LEN_PAT      = re.compile(r'len(\d+)', re.IGNORECASE)
PARENT_TAIL  = re.compile(r'(?:^|[_-])([MP])$')

def to_float(x: str):
    try: return float(x)
    except: return None

def is_pass(label: str) -> bool:
    return ("PASS" in (label or ""))

def is_lowdepth(label: str) -> bool:
    return ("low_depth" in (label or "").lower())

def is_uncallable(label: str) -> bool:
    return ("uncallable" in (label or "").lower())

def is_uncallable_e3(label: str) -> bool:
    return ("uncallable_stagee3" in (label or "").lower())

def read_table(path: str):
    with open(path, "r") as f:
        rdr = csv.reader(f, delimiter="\t", lineterminator="\n")
        headers = next(rdr)
        rows = [r for r in rdr]
    idx = {h:i for i,h in enumerate(headers)}
    return headers, rows, idx

def looks_like_header(row: List[str]) -> bool:
    if not row: return True
    j = "\t".join(c.lower() for c in row)
    return ("allele" in j) or ("contig" in j) or ("id" in j) or \
           (len(row)>=2 and (row[0].lower() in {"allele1","contig1","father"} or row[1].lower() in {"allele2","contig2","mother"}))

def read_pairs(path: str):
    pairs = []
    if not path: return pairs
    with open(path, "r") as f:
        rdr = csv.reader(f, delimiter="\t", lineterminator="\n")
        first = True
        for row in rdr:
            if not row: continue
            if first and looks_like_header(row):
                first = False; continue
            first = False
            if len(row)<2: continue
            a, b = row[0].strip(), row[1].strip()
            if a and b:
                pairs.append(tuple(sorted((a,b))))
    return pairs

def read_parental_homology(path: str):
    pairs = []
    with open(path, "r") as f:
        rdr = csv.reader(f, delimiter="\t", lineterminator="\n")
        first = True
        for row in rdr:
            if not row: continue
            if first and looks_like_header(row):
                first = False; continue
            first = False
            if len(row)<2: continue
            father, mother = row[0].strip(), row[1].strip()
            if father and mother:
                pairs.append((father, mother))
    return pairs

def parent_origin_from_id(cid: str) -> str:
    m = BUBBLE_PAT.match(cid)
    if m: return m.group(3)
    m2 = PARENT_TAIL.search(cid)
    if m2: return m2.group(1)
    return "U"

def main():
    ap = argparse.ArgumentParser(description="Callable genome with uncallableE3 regions and end-trimming")
    ap.add_argument("-i","--input", required=True)
    ap.add_argument("-o","--output-prefix", required=True)
    ap.add_argument("--parental-homology", required=True)
    ap.add_argument("-n","--bubblelike", default=None)
    ap.add_argument("--len-col", default=None)
    ap.add_argument("--cid-col", default="ContigID")
    ap.add_argument("--inh-col", default="inheritance")
    ap.add_argument("--mp-col",  default="MP_hybrid")
    ap.add_argument("--sim-thresh", type=float, default=0.5)
    ap.add_argument("--sim-op", choices=["gt","ge"], default="gt")
    ap.add_argument("--trim", type=int, default=100)
    args = ap.parse_args()

    headers, rows, idx = read_table(args.input)
    i_cid = idx.get(args.cid_col); i_inh = idx.get(args.inh_col); i_mp = idx.get(args.mp_col)
    if None in (i_cid, i_inh, i_mp):
        sys.exit("[ERROR] Missing required columns: ContigID / inheritance / MP_hybrid")

    i_len = idx.get(args.len_col) if args.len_col else None

    def full_len(r):
        if i_len is not None:
            try: return int(float(r[i_len]))
            except: return 0
        m = LEN_PAT.search(r[i_cid])
        if m: return int(m.group(1))
        return 0

    TRIM = max(0, int(args.trim))
    def trimmed_interval(L: int):
        s = TRIM; e = L - TRIM
        if e > s: return (s, e, e - s)
        return None

    id2row = { r[i_cid]:k for k,r in enumerate(rows) }

    def mp(cid):
        i = id2row.get(cid)
        try: return float(rows[i][i_mp]) if i is not None else None
        except: return None

    def inh_of(cid):
        i = id2row.get(cid)
        return None if i is None else rows[i][i_inh]

    def similar(mpv: float) -> bool:
        if mpv is None: return False
        return (mpv > args.sim_thresh) if args.sim_op=="gt" else (mpv >= args.sim_thresh)

    # parental homologous pairs
    parental_pairs = read_parental_homology(args.parental_homology)
    parental_pairs = [ (fa,mo) for (fa,mo) in parental_pairs if (fa in id2row and mo in id2row) ]

    # 1) excluded (both similar AND both inherited)
    excluded_pairs = []
    excluded_contigs = set()
    for fa, mo in parental_pairs:
        mfa, mmo = mp(fa), mp(mo)
        if mfa is None or mmo is None: continue
        if similar(mfa) and similar(mmo) and is_pass(inh_of(fa)) and is_pass(inh_of(mo)):
            excluded_pairs.append([fa, mo, inh_of(fa), inh_of(mo), mfa, mmo, "both_similar_and_inherited"])
            excluded_contigs.update([fa, mo])

    # 2) uncallable_stageE3 regions (both labels contain 'uncallable_stageE3')
    e3_pairs = []
    e3_contigs = set()
    for fa, mo in parental_pairs:
        ia, ib = inh_of(fa), inh_of(mo)
        if ia is None or ib is None: continue
        if is_uncallable_e3(ia) and is_uncallable_e3(ib):
            e3_pairs.append([fa, mo, ia, ib, mp(fa), mp(mo), "both_uncallable_stageE3"])
            e3_contigs.update([fa, mo])

    # 3) callable BED (post-trim)
    callable_bed = []
    for r in rows:
        cid = r[i_cid]; inh = r[i_inh]
        if is_pass(inh) and (not is_lowdepth(inh)) and (not is_uncallable(inh)) and (cid not in excluded_contigs):
            L = full_len(r)
            ti = trimmed_interval(L)
            if ti is not None:
                s,e,_ = ti
                callable_bed.append((cid, s, e))

    # 4) uncallableE3 BED (post-trim)
    e3_bed = []
    for cid in sorted(e3_contigs):
        r = rows[id2row[cid]]
        L = full_len(r)
        ti = trimmed_interval(L)
        if ti is not None:
            s,e,_ = ti
            e3_bed.append((cid, s, e))

    # write BEDs
    with open(args.output_prefix + ".callable_filtered.bed", "w", newline="") as f:
        for cid, s, e in callable_bed:
            f.write(f"{cid}\t{s}\t{e}\n")

    with open(args.output_prefix + ".uncallableE3.bed", "w", newline="") as f:
        for cid, s, e in e3_bed:
            f.write(f"{cid}\t{s}\t{e}\n")

    # pair details
    with open(args.output_prefix + ".callable_excluded_pairs.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["father_contig","mother_contig","father_inh","mother_inh","father_MP","mother_MP","reason"])
        for rec in excluded_pairs:
            w.writerow(rec)

    with open(args.output_prefix + ".uncallableE3_pairs.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["father_contig","mother_contig","father_inh","mother_inh","father_MP","mother_MP","reason"])
        for rec in e3_pairs:
            w.writerow(rec)

    # category stats (post-trim)
    def is_bubble_cid(cid: str) -> bool:
        return bool(BUBBLE_PAT.match(cid))

    bubblelike_pairs = read_pairs(args.bubblelike)
    bubblelike_nodes = set([x for p in bubblelike_pairs for x in p])

    mp_map  = { r[i_cid]: to_float(r[i_mp]) for r in rows }

    # full + trimmed length maps
    full_len_map = {}
    for r in rows:
        cid = r[i_cid]
        if i_len is not None:
            try:
                L = int(float(r[i_len]))
            except:
                m = LEN_PAT.search(cid); L = int(m.group(1)) if m else 0
        else:
            m = LEN_PAT.search(cid); L = int(m.group(1)) if m else 0
        full_len_map[cid] = L

    trimmed_len_map: Dict[str, int] = {}
    for cid, L in full_len_map.items():
        if L - args.trim > args.trim:
            trimmed_len_map[cid] = (L - args.trim) - args.trim
        else:
            trimmed_len_map[cid] = 0

    def parent_of(cid: str) -> str:
        m = BUBBLE_PAT.match(cid)
        if m: return m.group(3)
        m2 = PARENT_TAIL.search(cid)
        if m2: return m2.group(1)
        return "U"

    def similar_cid(cid: str) -> bool:
        mpv = mp_map.get(cid)
        if mpv is None: return False
        return (mpv > args.sim_thresh) if args.sim_op=="gt" else (mpv >= args.sim_thresh)

    cat_defs = [
        ("all",                      lambda cid: True),
        ("father_all",               lambda cid: parent_of(cid) == "P"),
        ("mother_all",               lambda cid: parent_of(cid) == "M"),
        ("unknown_parent",           lambda cid: parent_of(cid) == "U"),
        ("bubble_all",               lambda cid: is_bubble_cid(cid)),
        ("nonbubble_all",            lambda cid: not is_bubble_cid(cid)),
        ("primary",                  lambda cid: cid.startswith("primary_bubble")),
        ("secondary",                lambda cid: cid.startswith("secondary_bubble")),
        ("bubblelike",               lambda cid: (cid in bubblelike_nodes)),
        ("bubble_similar",           lambda cid: is_bubble_cid(cid) and similar_cid(cid)),
        ("nonbubble_similar",        lambda cid: (not is_bubble_cid(cid)) and similar_cid(cid)),
    ]

    excl_counts = defaultdict(lambda: [0,0])  # excluded (post-trim)
    call_counts = defaultdict(lambda: [0,0])  # callable (post-trim)
    e3_counts   = defaultdict(lambda: [0,0])  # uncallableE3 (post-trim)

    for cid in set(excluded_contigs):
        tl = trimmed_len_map.get(cid, 0)
        if tl <= 0: continue
        for cname, pred in cat_defs:
            if pred(cid):
                excl_counts[cname][0] += 1
                excl_counts[cname][1] += tl

    for cid, s, e in callable_bed:
        tl = e - s
        for cname, pred in cat_defs:
            if pred(cid):
                call_counts[cname][0] += 1
                call_counts[cname][1] += tl

    for cid, s, e in e3_bed:
        tl = e - s
        for cname, pred in cat_defs:
            if pred(cid):
                e3_counts[cname][0] += 1
                e3_counts[cname][1] += tl

    with open(args.output_prefix + ".callable_filter_stats.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["category","excluded_n","excluded_bp_posttrim","uncallableE3_n","uncallableE3_bp_posttrim","callable_n","callable_bp_posttrim","trim_bp_each_end"])
        for cname, _ in cat_defs:
            en, ebp = excl_counts[cname]
            e3n, e3bp = e3_counts[cname]
            cn, cbp = call_counts[cname]
            w.writerow([cname, en, ebp, e3n, e3bp, cn, cbp, args.trim])

if __name__ == "__main__":
    main()
