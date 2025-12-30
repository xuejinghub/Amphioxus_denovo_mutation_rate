#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
RBH on PAF BLOCKS only (no chaining), enforce uniqueness at BLOCK level only.
Aggregate passing blocks into a parental homologous contig list WITHOUT any contig-level uniqueness filtering.

Inputs:
  - paf_M2P: PAF from maternal (query) -> paternal (target)
  - paf_P2M: PAF from paternal (query) -> maternal (target)

Outputs:
  - <prefix>.block_pairs.tsv  : block-level reciprocal-best pairs (block-level uniqueness only)
  - <prefix>.contig_pairs.tsv : aggregated (M_contig, P_contig) stats without contig-level uniqueness filtering
"""

import sys, argparse
from collections import defaultdict, namedtuple

PAFBlock = namedtuple("PAFBlock",
    "q qlen qs qe strand t tlen ts te nmatch alnlen mapq")

def parse_paf(path, exclude_self=False):
    out = []
    with open(path) as f:
        for ln in f:
            if not ln.strip(): continue
            a = ln.rstrip("\n").split("\t")
            if len(a) < 12: continue
            q, qlen, qs, qe = a[0], int(a[1]), int(a[2]), int(a[3])
            strand = a[4]
            t, tlen, ts, te = a[5], int(a[6]), int(a[7]), int(a[8])
            nmatch, alnlen, mapq = int(a[9]), int(a[10]), int(a[11])
            if exclude_self and q == t:
                continue
            out.append(PAFBlock(q, qlen, qs, qe, strand, t, tlen, ts, te, nmatch, alnlen, mapq))
    return out

def merge_intervals(ivls):
    if not ivls: return [], 0
    ivls = sorted(ivls)
    out = []
    cs, ce = ivls[0]
    for s,e in ivls[1:]:
        if s <= ce:
            ce = max(ce, e)
        else:
            out.append((cs, ce)); cs, ce = s, e
    out.append((cs, ce))
    total = sum(e - s for s,e in out)
    return out, total

def overlap_len(a, b):
    s = max(a[0], b[0]); e = min(a[1], b[1])
    return max(0, e - s)

def overlap_frac(a, b):
    ov = overlap_len(a, b)
    la = a[1]-a[0]; lb = b[1]-b[0]
    den = min(la, lb) if min(la, lb) > 0 else 1
    return ov / den

# ---------- block-level 1->many guards (only at block level, not contig level) ----------

def block_overlaps_in_same_query(block, blocks_by_q, min_frac=0.5):
    """Does this block's query interval overlap (>=min_frac) another block mapping to a different target contig?"""
    q_rng = (block.qs, block.qe)
    arr = blocks_by_q.get(block.q, [])
    for g in arr:
        if g is block: continue
        if g.t == block.t: continue
        if overlap_frac(q_rng, (g.qs, g.qe)) >= min_frac:
            return True
    return False

def block_overlaps_in_same_target(block, blocks_by_t, min_frac=0.5):
    """Symmetric check on target axis."""
    t_rng = (block.ts, block.te)
    arr = blocks_by_t.get(block.t, [])
    for g in arr:
        if g is block: continue
        if g.q == block.q: continue
        if overlap_frac(t_rng, (g.ts, g.te)) >= min_frac:
            return True
    return False

def main():
    ap = argparse.ArgumentParser(description="Block-level RBH (no chaining), uniqueness only at block level; aggregate to contig list.")
    ap.add_argument("paf_M2P", help="PAF: maternal (query) -> paternal (target)")
    ap.add_argument("paf_P2M", help="PAF: paternal (query) -> maternal (target)")
    ap.add_argument("--out-prefix", default="rbh_blocks_simple", help="Output prefix")

    # block-level filters (no min-aln-block by request)
    ap.add_argument("--min-ident-block", type=float, default=0.95, help="Min identity for a single PAF block (matches/alnlen)")
    ap.add_argument("--overlap-frac", type=float, default=0.5, help="Min bidirectional overlap fraction for block-level RBH (on both axes)")
    ap.add_argument("--uni-overlap", type=float, default=0.5, help="Min same-axis overlap to flag 1->many ambiguity on a block")

    args = ap.parse_args()

    # 1) load PAF blocks (keep all)
    M2P_all = parse_paf(args.paf_M2P, exclude_self=False)
    P2M_all = parse_paf(args.paf_P2M, exclude_self=False)

    # 2) identity-only filter at block level (no min aligned length)
    def good_block(b):
        if b.alnlen <= 0: return False
        ident = b.nmatch / b.alnlen
        return ident >= args.min_ident_block

    M2P = [b for b in M2P_all if good_block(b)]
    P2M = [b for b in P2M_all if good_block(b)]

    # 3) indices for reverse lookup & ambiguity guards
    rev_index = defaultdict(list)   # (P_contig, M_contig) -> list[P2M blocks]
    for rb in P2M:
        rev_index[(rb.q, rb.t)].append(rb)

    by_q_M2P = defaultdict(list); by_t_M2P = defaultdict(list)
    for b in M2P:
        by_q_M2P[b.q].append(b); by_t_M2P[b.t].append(b)

    by_q_P2M = defaultdict(list); by_t_P2M = defaultdict(list)
    for b in P2M:
        by_q_P2M[b.q].append(b); by_t_P2M[b.t].append(b)

    # 4) forward: best reverse partner by bidirectional overlap; block-level ambiguity guards
    fwd_best = {}  # id(fb) -> (rb, score, fq, ft)
    for fb in M2P:
        if block_overlaps_in_same_query(fb, by_q_M2P, min_frac=args.uni_overlap):
            continue
        if block_overlaps_in_same_target(fb, by_t_M2P, min_frac=args.uni_overlap):
            continue
        cands = rev_index.get((fb.t, fb.q), [])
        if not cands: 
            continue
        fq_rng = (fb.qs, fb.qe)
        ft_rng = (fb.ts, fb.te)
        best = (None, 0.0, 0.0, 0.0)
        for rb in cands:
            rq_rng = (rb.qs, rb.qe)  # on paternal
            rt_rng = (rb.ts, rb.te)  # on maternal
            fq = overlap_frac(fq_rng, rt_rng)  # maternal axis
            ft = overlap_frac(ft_rng, rq_rng)  # paternal axis
            score = min(fq, ft)
            if score > best[1]:
                best = (rb, score, fq, ft)
        if best[1] >= args.overlap_frac:
            fwd_best[id(fb)] = best

    # 5) reverse: best forward partner for reciprocity; same ambiguity guards
    fwd_index = defaultdict(list)   # (M_contig, P_contig) -> list[M2P blocks]
    for fb in M2P:
        fwd_index[(fb.q, fb.t)].append(fb)

    rev_best = {}  # id(rb) -> (fb, score, fq, ft)
    for rb in P2M:
        if block_overlaps_in_same_query(rb, by_q_P2M, min_frac=args.uni_overlap):
            continue
        if block_overlaps_in_same_target(rb, by_t_P2M, min_frac=args.uni_overlap):
            continue
        cands = fwd_index.get((rb.t, rb.q), [])
        if not cands:
            continue
        rq_rng = (rb.qs, rb.qe)
        rt_rng = (rb.ts, rb.te)
        best = (None, 0.0, 0.0, 0.0)
        for fb in cands:
            fq = overlap_frac((fb.qs, fb.qe), rt_rng)
            ft = overlap_frac((fb.ts, fb.te), rq_rng)
            score = min(fq, ft)
            if score > best[1]:
                best = (fb, score, fq, ft)
        if best[1] >= args.overlap_frac:
            rev_best[id(rb)] = best

    # 6) reciprocal-best at block level only
    id_to_fb = {id(b): b for b in M2P}
    id_to_rb = {id(b): b for b in P2M}

    block_pairs = []  # (fb, rb, s_fw, s_rv)
    chosen_fwd_for_rb = { rid: rev_best[rid][0] for rid in rev_best }
    for fb in M2P:
        rec = fwd_best.get(id(fb))
        if not rec: continue
        rb, s_fw, fq_fw, ft_fw = rec
        rec_r = rev_best.get(id(rb))
        if not rec_r: continue
        if rec_r[0] is fb:
            block_pairs.append((fb, rb, s_fw, rec_r[1]))

    # 7) write block-level RBH table
    blk_out = args.out_prefix + ".block_pairs.tsv"
    with open(blk_out, "w") as fo:
        print("\t".join([
            "M_contig","M_len","M_qs","M_qe",
            "P_contig","P_len","P_ts","P_te",
            "strand_M2P","ident_M2P","aln_M2P","MAPQ_M2P",
            "ident_P2M","aln_P2M","MAPQ_P2M",
            "overlapScore_M2P","overlapScore_P2M"
        ]), file=fo)
        for fb, rb, s_fw, s_rv in block_pairs:
            print("\t".join(map(str,[
                fb.q, fb.qlen, fb.qs, fb.qe,
                rb.q, rb.qlen, rb.qs, rb.qe, 
                fb.strand,
                round(fb.nmatch/fb.alnlen,5), fb.alnlen, fb.mapq,
                round(rb.nmatch/rb.alnlen,5), rb.alnlen, rb.mapq,
                round(s_fw,3), round(s_rv,3)
            ])), file=fo)

    # 8) aggregate passing blocks -> (M_contig, P_contig) list (no contig-level uniqueness filtering)
    pair_stats = {}  # (M, P) -> aggregates
    M_ivls = defaultdict(list)
    P_ivls = defaultdict(list)

    for fb, rb, _, _ in block_pairs:
        M = fb.q; P = fb.t
        key = (M, P)
        st = pair_stats.get(key)
        if st is None:
            st = {"M_len": fb.qlen, "P_len": fb.tlen, "matches": 0, "aln": 0, "n_blocks": 0}
            pair_stats[key] = st
        st["matches"] += fb.nmatch
        st["aln"]     += fb.alnlen
        st["n_blocks"] += 1
        M_ivls[key].append((fb.qs, fb.qe))
        P_ivls[key].append((fb.ts, fb.te))

    # compute coverage & weighted identity
    for key, st in pair_stats.items():
        (M, P) = key
        _, Mcov = merge_intervals(M_ivls[key])
        _, Pcov = merge_intervals(P_ivls[key])
        st["M_cov_frac"] = Mcov / st["M_len"] if st["M_len"] > 0 else 0.0
        st["P_cov_frac"] = Pcov / st["P_len"] if st["P_len"] > 0 else 0.0
        st["ident_w"]    = st["matches"]/st["aln"] if st["aln"] > 0 else 0.0
        st["score"]      = st["matches"]

    outp = args.out_prefix + ".contig_pairs.tsv"
    with open(outp, "w") as fo:
        print("\t".join([
            "M_contig","M_len","P_contig","P_len",
            "identity_weighted","M_coverage","P_coverage",
            "sum_aligned","n_blocks"
        ]), file=fo)
        for (M,P), st in sorted(pair_stats.items(), key=lambda x: x[1]["score"], reverse=True):
            print("\t".join(map(str,[
                M, st["M_len"], P, st["P_len"],
                round(st["ident_w"],5),
                round(st["M_cov_frac"],4), round(st["P_cov_frac"],4),
                st["aln"], st["n_blocks"]
            ])), file=fo)

    sys.stderr.write(f"[done] block-level RBH pairs: {len(block_pairs)} ; contig pairs (aggregated, no uniqueness filtering): {len(pair_stats)}\n")

if __name__ == "__main__":
    main()
