#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, csv, re, argparse
from collections import defaultdict
from typing import Optional, Tuple, Set, List

BUBBLE_PAT = re.compile(r'^(primary|secondary)_bubble(\d+).*?([MP])$')
LEN_PAT    = re.compile(r'len(\d+)', re.IGNORECASE)

def is_pass(label: str) -> bool:
    return ("PASS" in (label or ""))

def is_lowdepth(label: str) -> bool:
    return ("low_depth" in (label or "").lower())

def to_float(x: str) -> Optional[float]:
    try:
        return float(x)
    except:
        return None

def read_table(path: str):
    with open(path, "r") as f:
        rdr = csv.reader(f, delimiter="\t")
        headers = next(rdr)
        rows = [r for r in rdr]
    return headers, rows, {h:i for i,h in enumerate(headers)}

def find_bubble_pairs(rows, i_cid) -> Set[Tuple[str,str]]:
    groups = defaultdict(dict)  # (num,src) -> {'primary': id, 'secondary': id}
    for row in rows:
        cid = row[i_cid]
        m = BUBBLE_PAT.match(cid)
        if m:
            role, num, src = m.groups()
            groups[(num,src)][role] = cid
    edges = set()
    for (_, _), members in groups.items():
        a = members.get("primary"); b = members.get("secondary")
        if a and b:
            edges.add(tuple(sorted((a,b))))
    return edges

def looks_like_header(row: List[str]) -> bool:
    if not row: return True
    j = "\t".join(c.lower() for c in row)
    return ("allele" in j) or ("contig" in j) or ("id" in j) or \
           (len(row)>=2 and (row[0].lower() in {"allele1","contig1"} or row[1].lower() in {"allele2","contig2"}))

def read_pairs(path: Optional[str]) -> List[Tuple[str,str]]:
    pairs = []
    if not path: return pairs
    with open(path, "r") as f:
        rdr = csv.reader(f, delimiter="\t")
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

def read_list(path: Optional[str]) -> Set[str]:
    s = set()
    if not path: return s
    with open(path, "r") as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.lower().startswith(("contig","allele","id")):
                continue
            s.add(line.split()[0])
    return s

def enforce_one_to_one(pairs: List[Tuple[str,str]]) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str]]]:
    count = defaultdict(int)
    for a,b in pairs:
        count[a]+=1; count[b]+=1
    kept, dropped = [], []
    for a,b in pairs:
        if count[a]==1 and count[b]==1:
            kept.append((a,b))
        else:
            dropped.append((a,b))
    return kept, dropped

def main():
    ap = argparse.ArgumentParser(description="StageA + pair-wise processing with revised Stage D & Stage E logic")
    ap.add_argument("-i","--input", required=True)
    ap.add_argument("-o","--output", required=True)
    ap.add_argument("-n","--nonbubble", dest="bubblelike", default=None, help="two-column TSV bubblelike 1:1 candidates")
    ap.add_argument("--bubblelike", dest="bubblelike2", default=None)
    ap.add_argument("-r","--real-nonbubble", dest="real_nonbubble", default=None, help="one-column TSV of REAL non-bubble contigs")
    ap.add_argument("--cid-col", default="ContigID")
    ap.add_argument("--inh-col", default="inheritance")
    ap.add_argument("--mp-col",  default="MP_hybrid")
    ap.add_argument("--len-col", default=None)
    ap.add_argument("--sim-thresh", type=float, default=0.5)
    ap.add_argument("--enforce-one-to-one", action="store_true", default=True)
    ap.add_argument("--orphan-label", default="nonInherited_stageF")
    ap.add_argument("--drop_two_similar", action="store_true", default=False,
                    help="Stage E low+low -> both uncallable_stageE3")
    args = ap.parse_args()

    bl_path = args.bubblelike or args.bubblelike2

    headers, rows, idx = read_table(args.input)
    i_cid = idx.get(args.cid_col); i_inh = idx.get(args.inh_col); i_mp = idx.get(args.mp_col)
    if None in (i_cid, i_inh, i_mp):
        sys.exit("[ERROR] Missing required columns: ContigID / inheritance / MP_hybrid")

    # length accessor
    i_len = None
    if args.len_col:
        i_len = idx.get(args.len_col)
        if i_len is None:
            sys.exit(f"[ERROR] --len-col '{args.len_col}' not found in header")

    def get_len(r):
        if i_len is not None:
            try: return int(float(r[i_len]))
            except: return 0
        m = LEN_PAT.search(r[i_cid])
        if m: return int(m.group(1))
        return 0

    kept_rows = rows[:]
    id2row = { r[i_cid]:k for k,r in enumerate(kept_rows) }

    # bubble & bubblelike pairs
    bubble_pairs     = find_bubble_pairs(kept_rows, i_cid)
    bubblelike_pairs_input = read_pairs(bl_path)
    bubblelike_pairs_input = [p for p in bubblelike_pairs_input if (p[0] in id2row and p[1] in id2row)]

    dropped_pairs = []
    if args.enforce_one_to_one and bubblelike_pairs_input:
        bubblelike_pairs, dropped_pairs = enforce_one_to_one(bubblelike_pairs_input)
    else:
        bubblelike_pairs = bubblelike_pairs_input

    all_pairs = set(bubble_pairs) | set(bubblelike_pairs)
    pair_type = {}
    for a,b in all_pairs:
        t = []
        if (a,b) in bubble_pairs:     t.append("bubble")
        if (a,b) in bubblelike_pairs: t.append("bubblelike")
        pair_type[(a,b)] = ",".join(t) if t else "unknown"

    # Stage stats tracker: [inh_n, inh_bp, non_n, non_bp, unc_n, unc_bp]
    stage_stats = defaultdict(lambda: [0,0,0,0,0,0])
    def add_stage(stage, label_assigned, length):
        l = (label_assigned or "").lower()
        if "uncallable" in l:
            stage_stats[stage][4] += 1
            stage_stats[stage][5] += length
        elif ("pass" in (label_assigned or "")) and ("uncallable" not in l):
            stage_stats[stage][0] += 1
            stage_stats[stage][1] += length
        else:
            stage_stats[stage][2] += 1
            stage_stats[stage][3] += length

    special_rows = []
    uncallable_rows = []

    # ---------------- Stage A: real-nonbubble ----------------
    real_set = read_list(args.real_nonbubble)
    for cid in sorted(real_set):
        ri = id2row.get(cid)
        if ri is None: continue
        old = kept_rows[ri][i_inh]
        ln  = get_len(kept_rows[ri])
        if is_lowdepth(old):
            special_rows.append(["stageA_real_nonbubble_keep_low_depth", "stageA", cid, old, kept_rows[ri][i_mp], "-", "-", old, "keep"])
            add_stage("A", old, ln)
        elif is_pass(old):
            kept_rows[ri][i_inh] = "PASS_stageA"
            special_rows.append(["stageA_real_nonbubble_keep_PASS", "stageA", cid, old, kept_rows[ri][i_mp], "-", "-", old, "keep"])
            add_stage("A", old, ln)
        else:
            kept_rows[ri][i_inh] = "PASS_stageA"
            special_rows.append(["stageA_real_nonbubble_set_PASS", "stageA", cid, old, kept_rows[ri][i_mp], "-", "-", kept_rows[ri][i_inh], "set_PASS_stageA"])
            add_stage("A", "PASS_stageA", ln)

    # ---------------- Stage B: both low_depth -> uncallable (log only) ----------------
    for a,b in sorted(all_pairs):
        ia = id2row.get(a); ib = id2row.get(b)
        if ia is None or ib is None: continue
        inh_a = kept_rows[ia][i_inh]; inh_b = kept_rows[ib][i_inh]
        if is_lowdepth(inh_a) and is_lowdepth(inh_b):
            uncallable_rows.append(["stageB_both_low_depth_uncallable", pair_type[(a,b)], a, kept_rows[ia][i_mp], b, kept_rows[ib][i_mp]])
            add_stage("B", "uncallable", get_len(kept_rows[ia])); add_stage("B", "uncallable", get_len(kept_rows[ib]))
            special_rows.append(["stageB_both_low_depth_keep", pair_type[(a,b)], a, inh_a, kept_rows[ia][i_mp], b, inh_b, kept_rows[ib][i_mp], "log_uncallable"])

    # ---------------- Stage C: single low_depth -> mate PASS_stageC ----------------
    for a,b in sorted(all_pairs):
        ia = id2row.get(a); ib = id2row.get(b)
        if ia is None or ib is None: continue
        inh_a = kept_rows[ia][i_inh]; inh_b = kept_rows[ib][i_inh]
        ld_a, ld_b = is_lowdepth(inh_a), is_lowdepth(inh_b)
        if ld_a ^ ld_b:
            if ld_a and not ld_b:
                mate = b; mi = ib; ld_cid = a; li = ia
            else:
                mate = a; mi = ia; ld_cid = b; li = ib
            old_mate = kept_rows[mi][i_inh]
            kept_rows[mi][i_inh] = "PASS_stageC"
            special_rows.append(["stageC_single_low_depth_force_PASS_stageC", pair_type[(a,b)], mate, old_mate, kept_rows[mi][i_mp], ld_cid, kept_rows[li][i_inh], "PASS_stageC", "assign_PASS_stageC"])
            add_stage("C", "PASS_stageC", get_len(kept_rows[mi])); add_stage("C", kept_rows[li][i_inh], get_len(kept_rows[li]))

    def mp(cid):
        i = id2row.get(cid)
        return None if i is None else to_float(kept_rows[i][i_mp])

    # ---------------- Stage D: revised high-priority handling ----------------
    for a,b in sorted(all_pairs):
        ia = id2row.get(a); ib = id2row.get(b)
        if ia is None or ib is None: continue
        inh_a = kept_rows[ia][i_inh]; inh_b = kept_rows[ib][i_inh]
        mpa, mpb = mp(a), mp(b)

        if is_lowdepth(inh_a) or is_lowdepth(inh_b):
            continue

        high_a = is_pass(inh_a) and (mpa is not None and mpa < args.sim_thresh)
        high_b = is_pass(inh_b) and (mpb is not None and mpb < args.sim_thresh)

        # exactly one high-priority
        if high_a ^ high_b:
            if high_a:
                hp, op = a, b
            else:
                hp, op = b, a
            hi = id2row[hp]; oi = id2row[op]
            old_h, old_o = kept_rows[hi][i_inh], kept_rows[oi][i_inh]
            kept_rows[hi][i_inh] = "PASS_stageD"
            kept_rows[oi][i_inh] = "nonInherited_stageD"
            special_rows.append(["stageD_single_high_priority", pair_type[(a,b)], hp, old_h, kept_rows[hi][i_mp], op, old_o, kept_rows[oi][i_inh], "hp->PASS_stageD, other->nonInherited_stageD"])
            add_stage("D", "PASS_stageD", get_len(kept_rows[hi])); add_stage("D", "nonInherited_stageD", get_len(kept_rows[oi]))
            continue

        # both high-priority
        if high_a and high_b:
            if (mpa < mpb) or (mpa == mpb and a < b):
                ctrl, other = a, b
            else:
                ctrl, other = b, a
            ci = id2row[ctrl]; oi = id2row[other]
            old_c, old_o = kept_rows[ci][i_inh], kept_rows[oi][i_inh]
            kept_rows[ci][i_inh] = "PASS_stageD"
            kept_rows[oi][i_inh] = "nonInherited_stageD"
            special_rows.append(["stageD_both_high_priority", pair_type[(a,b)], ctrl, old_c, kept_rows[ci][i_mp], other, old_o, kept_rows[oi][i_inh], "lowerMP->PASS_stageD"])
            add_stage("D", "PASS_stageD", get_len(kept_rows[ci])); add_stage("D", "nonInherited_stageD", get_len(kept_rows[oi]))
            continue

    # ---------------- Stage E: mid/low combos (updated logic) ----------------
    for a,b in sorted(all_pairs):
        ia = id2row.get(a); ib = id2row.get(b)
        if ia is None or ib is None: continue

        inh_a = kept_rows[ia][i_inh]; inh_b = kept_rows[ib][i_inh]
        mpa, mpb = mp(a), mp(b)

        # skip if any high-priority PASS already handled by D
        high_a = is_pass(inh_a) and (mpa is not None and mpa < args.sim_thresh)
        high_b = is_pass(inh_b) and (mpb is not None and mpb < args.sim_thresh)
        if high_a or high_b:
            continue
        if is_lowdepth(inh_a) or is_lowdepth(inh_b):
            continue

        mid_a = (mpa is not None and mpa < args.sim_thresh)
        mid_b = (mpb is not None and mpb < args.sim_thresh)
        low_a = (mpa is not None and mpa >= args.sim_thresh)
        low_b = (mpb is not None and mpb >= args.sim_thresh)

        # E1: mid+mid
        if mid_a and mid_b:
            if (mpa < mpb) or (mpa == mpb and a < b):
                noni, passi = a, b
            else:
                noni, passi = b, a
            ni = id2row[noni]; pi = id2row[passi]
            old_n = kept_rows[ni][i_inh]; old_p = kept_rows[pi][i_inh]
            kept_rows[ni][i_inh] = "nonInherited"
            kept_rows[pi][i_inh] = "PASS_stageE1"
            special_rows.append(["stageE_mid_mid_split", pair_type[(a,b)], noni, old_n, kept_rows[ni][i_mp], passi, old_p, kept_rows[pi][i_inh], "nonInherited_vs_PASS_stageE1"])
            add_stage("E", "nonInherited", get_len(kept_rows[ni])); add_stage("E", "PASS_stageE1", get_len(kept_rows[pi]))
            continue

        # E2: mid+low
        if (mid_a and low_b) or (low_a and mid_b):
            if mid_a and low_b:
                mid_c, low_c = a, b
            else:
                mid_c, low_c = b, a
            mi = id2row[mid_c]; li2 = id2row[low_c]
            old_m = kept_rows[mi][i_inh]; old_l = kept_rows[li2][i_inh]
            kept_rows[mi][i_inh]  = "nonInherited_stageE2"
            kept_rows[li2][i_inh] = "PASS_stageE2"
            special_rows.append(["stageE_mid_low_split", pair_type[(a,b)], mid_c, old_m, kept_rows[mi][i_mp], low_c, old_l, kept_rows[li2][i_inh], "mid_nonInherited_stageE2_low_PASS_stageE2"])
            add_stage("E", "nonInherited_stageE2", get_len(kept_rows[mi])); add_stage("E", "PASS_stageE2", get_len(kept_rows[li2]))
            continue

        # E3: low+low
        if low_a and low_b:
            old_a, old_b = kept_rows[ia][i_inh], kept_rows[ib][i_inh]
            if args.drop_two_similar:
                kept_rows[ia][i_inh] = "uncallable_stageE3"
                kept_rows[ib][i_inh] = "uncallable_stageE3"
                special_rows.append(["stageE_low_low_drop_two_similar", pair_type[(a,b)], a, old_a, kept_rows[ia][i_mp], b, old_b, kept_rows[ib][i_inh], "both_uncallable_stageE3"])
                add_stage("E", "uncallable_stageE3", get_len(kept_rows[ia])); add_stage("E", "uncallable_stageE3", get_len(kept_rows[ib]))
            else:
                # controller = lower MP (tie -> lexicographic), then use its CURRENT PASS status
                if (mpa < mpb) or (mpa == mpb and a < b):
                    controller, other = a, b
                else:
                    controller, other = b, a
                ci = id2row[controller]; oi = id2row[other]
                old_c = kept_rows[ci][i_inh]; old_o = kept_rows[oi][i_inh]
                if is_pass(old_c):
                    kept_rows[ci][i_inh] = "PASS_stageE3"
                    kept_rows[oi][i_inh] = "nonInherited_stageE3"
                    reason = "low+low_use_lowerMP_currentPASS"
                else:
                    kept_rows[ci][i_inh] = "nonInherited_stageE3"
                    kept_rows[oi][i_inh] = "PASS_stageE3"
                    reason = "low+low_use_lowerMP_currentNonPASS"
                special_rows.append(["stageE_low_low_controller_by_current_label", pair_type[(a,b)], controller, old_c, kept_rows[ci][i_mp], other, old_o, kept_rows[oi][i_inh], reason])
                add_stage("E", kept_rows[ci][i_inh], get_len(kept_rows[ci])); add_stage("E", kept_rows[oi][i_inh], get_len(kept_rows[oi]))

    # ---------------- Stage F: non-bubble orphan similar -> not inherited ----------------
    nb_kept_nodes = set()
    for a,b in bubblelike_pairs:
        nb_kept_nodes.add(a); nb_kept_nodes.add(b)

    def is_bubble_cid(cid: str) -> bool:
        return bool(BUBBLE_PAT.match(cid))

    for cid, ri in list(id2row.items()):
        if cid in nb_kept_nodes: continue
        if is_bubble_cid(cid): continue
        if cid in real_set:    continue
        inh = kept_rows[ri][i_inh]
        mpv = to_float(kept_rows[ri][i_mp])
        if mpv is None: continue
        if mpv >= args.sim_thresh:
            if (not is_lowdepth(inh)) and (inh != "unsolved_comfilct"):
                before = kept_rows[ri][i_inh]
                kept_rows[ri][i_inh] = args.orphan_label
                special_rows.append(["stageF_orphan_similar_not_inherited", "bubblelike_orphan",
                                     cid, before, kept_rows[ri][i_mp], "-", "-", kept_rows[ri][i_inh], "set_orphan_not_inherited"])
                add_stage("F", args.orphan_label, get_len(kept_rows[ri]))

    # ---------------- Outputs ----------------
    out_path = args.output
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(headers); w.writerows(kept_rows)

    with open(out_path + ".special_cases.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["scenario","pair_type_or_stage","source","source_inh_before","source_MP","target","target_inh_before","target_inh_after","action"])
        for rec in special_rows: w.writerow(rec)

    with open(out_path + ".uncallable_pairs.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["reason","pair_type","allele1","allele1_MP","allele2","allele2_MP"])
        for rec in uncallable_rows: w.writerow(rec)

    if dropped_pairs:
        with open(out_path + ".bubblelike_non_one_to_one.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["allele1","allele2","note"])
            for a,b in dropped_pairs:
                w.writerow([a,b,"dropped_non_one_to_one"])

    # ---------------- Final label × category stats ----------------
    bubblelike_nodes = nb_kept_nodes
    categories = [
        ("all",                         lambda cid, mp: True),
        ("bubble_all",                  lambda cid, mp: is_bubble_cid(cid)),
        ("nonbubble_excl_blike",       lambda cid, mp: (not is_bubble_cid(cid)) and (cid not in bubblelike_nodes)),
        ("primary",                    lambda cid, mp: cid.startswith("primary_bubble")),
        ("secondary",                  lambda cid, mp: cid.startswith("secondary_bubble")),
        ("bubblelike",                 lambda cid, mp: (cid in bubblelike_nodes)),
        ("bubble_similar",             lambda cid, mp: is_bubble_cid(cid) and (mp is not None and mp >= args.sim_thresh)),
        ("nonbubble_similar",          lambda cid, mp: (not is_bubble_cid(cid)) and (cid not in bubblelike_nodes) and (mp is not None and mp >= args.sim_thresh)),
    ]

    def get_len2(r):
        if i_len is not None:
            try: return int(float(r[i_len]))
            except: return 0
        m = LEN_PAT.search(r[i_cid])
        if m: return int(m.group(1))
        return 0

    agg = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    for r in kept_rows:
        cid = r[i_cid]
        inh = r[i_inh]
        mpv = to_float(r[i_mp])
        ln  = get_len2(r)
        for cname, pred in categories:
            if pred(cid, mpv):
                agg[inh][cname][0] += 1
                agg[inh][cname][1] += ln

    stats_path = out_path + ".stats_by_label.tsv"
    with open(stats_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        hdr = ["inheritance_label"]
        for cname, _ in categories:
            hdr += [f"{cname}_n", f"{cname}_bp"]
        w.writerow(hdr)
        for label in sorted(agg.keys()):
            row = [label]
            for cname, _ in categories:
                n, bp = agg[label][cname]
                row += [str(n), str(bp)]
            w.writerow(row)

    # ---------------- Stage summary stats ----------------
    stage_stats_path = out_path + ".stage_stats.tsv"
    with open(stage_stats_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["stage","inherited_n","inherited_bp","non_inherited_n","non_inherited_bp","uncallable_n","uncallable_bp"])
        for stg in ["A","B","C","D","E","F"]:
            n = stage_stats.get(stg, [0,0,0,0,0,0])
            w.writerow([stg] + [str(x) for x in n])

if __name__ == "__main__":
    main()
