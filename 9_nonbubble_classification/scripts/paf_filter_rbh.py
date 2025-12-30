#!/usr/bin/env python3
import sys, gzip, argparse

def open_any(path):
    return sys.stdin if path == "-" else (gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r"))

def merge_len(intervals):
    """intervals: list[(start,end)] 0-based half-open; return union length"""
    if not intervals:
        return 0
    intervals.sort()
    s, e = intervals[0]
    total = 0
    for a, b in intervals[1:]:
        if a > e:
            total += (e - s)
            s, e = a, b
        else:
            if b > e:
                e = b
    total += (e - s)
    return total

def overlaps(a1, a2, b1, b2):
    """half-open [a1,a2) vs [b1,b2) overlap"""
    return (a1 < b2) and (b1 < a2)

def parse_args():
    ap = argparse.ArgumentParser(
        description="Aggregate PAF by (q,t), filter by length-weighted identity and coverage; "
                    "coverage can be relaxed if any fragment touches either contig end within a window; "
                    "optionally compute fragment-level RBH with endpoint tolerance."
    )
    ap.add_argument("--paf", required=True, help="Input PAF (can be .gz). Use '-' for stdin.")
    ap.add_argument("--min-id", type=float, default=0.90, help="Minimum length-weighted identity (exclusive, default 0.90).")
    ap.add_argument("--min-cov", type=float, default=0.70, help="Minimum coverage (inclusive, default 0.70).")
    ap.add_argument("--end-relax", type=int, default=100,
                    help="If any fragment overlaps either end-window (both sides) within this many bp, "
                         "coverage filter is relaxed (default 100).")
    ap.add_argument("--out-pairs", default="filtered_pairs.tsv", help="Output TSV of filtered (q,t) pairs.")
    ap.add_argument("--out-failed-filter", default="failed_pairs.tsv", help="Output TSV of pairs failing the filter, with reasons.")
    ap.add_argument("--out-rbh", default="rbh.tsv", help="Output TSV of fragment-level RBH pairs (if RBH enabled).")
    ap.add_argument("--out-failed-rbh", default="failed_rbh.tsv", help="Output TSV of fragments that passed filters but failed RBH (if enabled).")
    ap.add_argument("--skip-rbh", action="store_true", help="Skip RBH computation/outputs.")
    ap.add_argument("--rbh-window", type=int, default=100,
                    help="Endpoint tolerance (bp) for fragment-level RBH (default 0).")
    return ap.parse_args()

def score_for_rbh_frag(hit):
    """Fragment-level ranking: prefer larger nmatch, then alnlen."""
    return (hit["nmatch"], hit["alnlen"])

def main():
    args = parse_args()
    end_relax = max(0, int(args.end_relax))

    # -------- Aggregate by (q,t) --------
    pairs = {}  # (q,t) -> dict: qlen,tlen, q_ints, t_ints, nmatch, alnlen, hits[]
    with open_any(args.paf) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q, qlen, qs, qe, strand, t, tlen, ts, te, nmatch, alnlen = (
                f[0], int(f[1]), int(f[2]), int(f[3]), f[4],
                f[5], int(f[6]), int(f[7]), int(f[8]), int(f[9]), int(f[10])
            )
            key = (q, t)
            rec = pairs.get(key)
            if rec is None:
                rec = {
                    "qlen": qlen, "tlen": tlen,
                    "q_ints": [], "t_ints": [],
                    "nmatch": 0, "alnlen": 0,
                    "hits": []
                }
                pairs[key] = rec
            else:
                if qlen > rec["qlen"]:
                    rec["qlen"] = qlen
                if tlen > rec["tlen"]:
                    rec["tlen"] = tlen
            rec["q_ints"].append((qs, qe))
            rec["t_ints"].append((ts, te))
            rec["nmatch"] += nmatch
            rec["alnlen"] += alnlen
            rec["hits"].append({
                "q": q, "qlen": qlen, "qs": qs, "qe": qe,
                "t": t, "tlen": tlen, "ts": ts, "te": te,
                "strand": strand,
                "nmatch": nmatch, "alnlen": alnlen
            })

    # -------- Compute stats & filter (with end-relax on coverage) --------
    filtered = {}   # (q,t) -> stats dict + hits
    with open(args.out_pairs, "w") as out_pass, open(args.out_failed_filter, "w") as out_fail:
        header = ["q","qlen","t","tlen","ID","coverage","union_q","union_t","union_min",
                  "sum_nmatch","sum_alnlen","cov_relaxed"]
        out_pass.write("\t".join(header) + "\n")
        out_fail.write("\t".join(header + ["fail_reason"]) + "\n")

        for (q,t), ag in pairs.items():
            qlen = ag["qlen"]; tlen = ag["tlen"]
            union_q = merge_len(ag["q_ints"])
            union_t = merge_len(ag["t_ints"])
            union_min = min(union_q, union_t)
            denom = min(qlen, tlen) if min(qlen, tlen) > 0 else 1
            coverage = union_min / denom
            ID = (ag["nmatch"] / ag["alnlen"]) if ag["alnlen"] > 0 else 0.0

            # 是否触发“端点放宽”
            cov_relaxed = 0
            if end_relax > 0:
                # 端窗口（clip 到 [0,len]）
                qL1, qL2 = 0, min(end_relax, qlen)
                qR1, qR2 = max(qlen - end_relax, 0), qlen
                tL1, tL2 = 0, min(end_relax, tlen)
                tR1, tR2 = max(tlen - end_relax, 0), tlen
                for h in ag["hits"]:
                    if (overlaps(h["qs"], h["qe"], qL1, qL2) or
                        overlaps(h["qs"], h["qe"], qR1, qR2) or
                        overlaps(h["ts"], h["te"], tL1, tL2) or
                        overlaps(h["ts"], h["te"], tR1, tR2)):
                        cov_relaxed = 1
                        break

            passes_id = (ID > args.min_id)
            passes_cov = (coverage >= args.min_cov) or (cov_relaxed == 1)

            base_cols = [
                q, str(qlen), t, str(tlen),
                f"{ID:.6f}", f"{coverage:.6f}",
                str(union_q), str(union_t), str(union_min),
                str(ag["nmatch"]), str(ag["alnlen"]), str(cov_relaxed)
            ]

            if passes_id and passes_cov:
                rec = {
                    "q": q, "qlen": qlen, "t": t, "tlen": tlen,
                    "ID": ID, "coverage": coverage,
                    "union_q": union_q, "union_t": union_t, "union_min": union_min,
                    "sum_nmatch": ag["nmatch"], "sum_alnlen": ag["alnlen"],
                    "cov_relaxed": cov_relaxed,
                    "hits": ag["hits"]
                }
                filtered[(q,t)] = rec
                out_pass.write("\t".join(base_cols) + "\n")
            else:
                reasons = []
                if ag["alnlen"] <= 0:
                    reasons.append("no_alnlen")
                if not passes_id:
                    reasons.append("id_below_threshold")
                if (cov_relaxed == 0) and (coverage < args.min_cov):
                    reasons.append("cov_below_threshold")
                out_fail.write("\t".join(base_cols + [";".join(reasons) if reasons else "unknown"]) + "\n")

    # -------- RBH (fragment-level, optional) --------
    if args.skip_rbh:
        print("[INFO] Skipping RBH (--skip-rbh).", file=sys.stderr)
        return

    hits_by_pair = {k: v["hits"] for k, v in filtered.items()}
    window = max(0, int(args.rbh_window))

    def match_within_window(h, r):
        return (abs(h["qs"] - r["ts"]) <= window and
                abs(h["qe"] - r["te"]) <= window and
                abs(h["ts"] - r["qs"]) <= window and
                abs(h["te"] - r["qe"]) <= window)

    rbh_pairs = []
    used_rev = set()
    failed_rbh = []

    for (q,t), rec in filtered.items():
        forward_hits = hits_by_pair.get((q,t), [])
        reverse_hits = hits_by_pair.get((t,q), [])
        if not reverse_hits:
            for h in forward_hits:
                failed_rbh.append({
                    "h": h, "reason": "no_reverse_pair_passed_filter",
                    "q_best_t": t, "t_best_q": "NA"
                })
            continue
        for h in forward_hits:
            best_r = None
            best_sc = None
            for r in reverse_hits:
                if match_within_window(h, r):
                    sc = score_for_rbh_frag(r)
                    if (best_r is None) or (sc > best_sc):
                        best_r, best_sc = r, sc
            if best_r is not None:
                keyr = (best_r["q"], best_r["qs"], best_r["qe"], best_r["t"], best_r["ts"], best_r["te"])
                if keyr in used_rev:
                    failed_rbh.append({"h": h, "reason": "reverse_fragment_already_paired",
                                       "q_best_t": t, "t_best_q": q})
                else:
                    used_rev.add(keyr)
                    rbh_pairs.append((h, best_r))
            else:
                failed_rbh.append({"h": h, "reason": "no_reciprocal_within_window",
                                   "q_best_t": t, "t_best_q": q})

    with open(args.out_rbh, "w") as out_rbh:
        header = ["q","qs","qe","qlen","t","ts","te","tlen","strand",
                  "nmatch","alnlen","ID_fragment","rbh_window","cov_relax_pair"]
        out_rbh.write("\t".join(header) + "\n")
        for h, r in rbh_pairs:
            idf = (h["nmatch"]/h["alnlen"]) if h["alnlen"]>0 else 0.0
            cov_relax_pair = filtered[(h["q"], h["t"])]["cov_relaxed"]
            out_rbh.write("\t".join([
                h["q"], str(h["qs"]), str(h["qe"]), str(h["qlen"]),
                h["t"], str(h["ts"]), str(h["te"]), str(h["tlen"]),
                h["strand"],
                str(h["nmatch"]), str(h["alnlen"]), f"{idf:.6f}",
                str(window), str(cov_relax_pair)
            ]) + "\n")

    with open(args.out_failed_rbh, "w") as out_fr:
        header = ["q","qs","qe","qlen","t","ts","te","tlen","strand",
                  "nmatch","alnlen","ID_fragment","fail_reason","rbh_window",
                  "q_best_t","t_best_q","cov_relax_pair"]
        out_fr.write("\t".join(header) + "\n")
        for item in failed_rbh:
            h = item["h"]
            idf = (h["nmatch"]/h["alnlen"]) if h["alnlen"]>0 else 0.0
            cov_relax_pair = filtered.get((h["q"], h["t"]), {}).get("cov_relaxed", 0)
            out_fr.write("\t".join([
                h["q"], str(h["qs"]), str(h["qe"]), str(h["qlen"]),
                h["t"], str(h["ts"]), str(h["te"]), str(h["tlen"]),
                h["strand"],
                str(h["nmatch"]), str(h["alnlen"]), f"{idf:.6f}",
                item.get("reason","not_reciprocal_best"), str(window),
                str(item.get("q_best_t","NA")), str(item.get("t_best_q","NA")),
                str(cov_relax_pair)
            ]) + "\n")

if __name__ == "__main__":
    main()
