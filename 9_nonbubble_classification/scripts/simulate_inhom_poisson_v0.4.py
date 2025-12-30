#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, argparse, random, math, gzip, bisect
from collections import defaultdict, namedtuple, Counter

try:
    import numpy as np
except ImportError:
    np = None

try:
    from scipy.stats import pearsonr
except ImportError:
    pearsonr = None

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None

def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    name, buf, order, seqs = None, [], [], {}
    with opener(path, "rt") as fh:
        for line in fh:
            if not line: continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                    buf = []
                name = line[1:].strip().split()[0]
                order.append(name)
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf)
    return order, seqs

def is_callable_base(b):
    return b in "ACGTacgt"

def count_callable(seq, start, end):
    s = seq[start:end]
    cnt = 0
    for ch in s:
        if ch in "ACGTacgt":
            cnt += 1
    return cnt

Window = namedtuple("Window", ["chrom","start","end","callable"])

def make_windows(seqs, order, win):
    per_chrom = {}
    all_windows = []
    for chrom in order:
        seq = seqs[chrom]
        L = len(seq)
        lst = []
        i = 0
        while i < L:
            j = min(L, i+win)
            c = count_callable(seq, i, j)
            w = Window(chrom, i, j, c)
            lst.append(w)
            all_windows.append(w)
            i = j
        per_chrom[chrom] = lst
    return per_chrom, all_windows

def load_observed_from_vcf(vcf_path, per_chrom_windows, win):
    opener = gzip.open if vcf_path.endswith(".gz") else open
    obs = {}
    for chrom, wlist in per_chrom_windows.items():
        obs[chrom] = [0]*len(wlist)
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if (not line) or line[0] == "#":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, pos, _id, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            if chrom not in per_chrom_windows:
                continue
            if len(ref) != 1:
                continue
            alt1 = alt.split(",")[0]
            if len(alt1) != 1:
                continue
            try:
                p = int(pos) - 1
            except:
                continue
            wlist = per_chrom_windows[chrom]
            if p < 0 or p >= wlist[-1].end:
                continue
            idx = p // win
            if 0 <= idx < len(wlist):
                obs[chrom][idx] += 1
    return obs

def load_observed_from_tsv(tsv_path, per_chrom_windows, win):
    obs = {chrom: [0]*len(wlist) for chrom, wlist in per_chrom_windows.items()}
    opener = gzip.open if tsv_path.endswith(".gz") else open
    with opener(tsv_path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            chrom, s, e, cnt = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
            if chrom not in per_chrom_windows:
                continue
            idx = s // win
            if 0 <= idx < len(per_chrom_windows[chrom]):
                obs[chrom][idx] += cnt
    return obs

def multinomial_allocate(total, probs):
    if total == 0:
        return [0]*len(probs)
    if np is not None:
        return list(np.random.multinomial(total, probs))
    raw = [p*total for p in probs]
    ints = [int(math.floor(x)) for x in raw]
    rem = total - sum(ints)
    fr = sorted(((raw[i]-ints[i], i) for i in range(len(raw))), reverse=True)
    for k in range(rem):
        ints[fr[k][1]] += 1
    return ints

def poisson_allocate(lams, rng):
    if np is not None:
        return list(np.random.poisson(lam=lams))
    out=[]
    for lam in lams:
        if lam<=0:
            out.append(0); continue
        L=math.exp(-lam); k=0; p=1.0
        while True:
            k+=1; p*=rng.random()
            if p<=L:
                out.append(k-1); break
    return out

def merge_intervals(intervals):
    if not intervals:
        return []
    intervals.sort()
    merged = [intervals[0]]
    for s,e in intervals[1:]:
        ls, le = merged[-1]
        if s <= le:
            merged[-1] = (ls, max(le, e))
        else:
            merged.append((s,e))
    return merged

def place_nonvariant_regions(seqs, order, total_count, rng, length_choices, max_tries=2000):
    """
    随机在各染色体放置 total_count 个非变异区（不重叠）。
    长度从 length_choices 中等概率抽取。
    返回: dict chrom -> list[(start,end,length)]
    """
    if total_count <= 0:
        return {chrom: [] for chrom in order}

    chrom_lengths = [len(seqs[c]) for c in order]
    genome_len = float(sum(chrom_lengths))
    probs = [L/genome_len for L in chrom_lengths]
    counts = multinomial_allocate(total_count, probs)

    out = {chrom: [] for chrom in order}
    for chrom, want in zip(order, counts):
        Lc = len(seqs[chrom])
        if Lc == 0 or want == 0:
            continue
        placed = []
        tries = 0
        while len(placed) < want and tries < max_tries:
            tries += 1
            length = rng.choice(length_choices)
            if length >= Lc:
                continue
            start = rng.randint(0, Lc - length)
            end = start + length
            overlap = False
            for s,e,_ in placed:
                if not (end <= s or start >= e):
                    overlap = True
                    break
            if overlap:
                continue
            placed.append((start,end,length))
        placed_intervals = merge_intervals([(s,e) for s,e,_ in placed])
        out[chrom] = [(s,e,e-s) for (s,e) in placed_intervals]
    return out

def build_interval_index(mask_intervals_by_chrom):
    idx = {}
    for chrom, ivals in mask_intervals_by_chrom.items():
        ivals_sorted = sorted([(s,e) for s,e,_L in ivals])
        starts = [s for s,_ in ivals_sorted]
        idx[chrom] = (ivals_sorted, starts)
    return idx

def point_in_intervals(intervals, starts, pos):
    if not intervals:
        return False
    i = bisect.bisect_right(starts, pos) - 1
    if i >= 0:
        s,e = intervals[i]
        if s <= pos < e:
            return True
    return False

def choose_positions_excluding(seq, start, end, k, rng, mask_index_for_chrom):
    intervals, starts = mask_index_for_chrom
    cand = []
    local = seq[start:end]
    for i, ch in enumerate(local):
        if not is_callable_base(ch):
            continue
        pos = start + i
        if not point_in_intervals(intervals, starts, pos):
            cand.append(pos)
    if not cand:
        return []
    if k >= len(cand):
        return cand[:]
    rng.shuffle(cand)
    return cand[:k]

def choose_positions_in_interval(seq, start, end, k, rng):
    """在 [start,end) 内从 A/C/G/T 位点中均匀选 k 个不重复位置。"""
    cand = []
    local = seq[start:end]
    for i, ch in enumerate(local):
        if is_callable_base(ch):
            cand.append(start + i)
    if not cand:
        return []
    if k >= len(cand):
        return cand[:]
    rng.shuffle(cand)
    return cand[:k]

def random_alt_base(ref, rng):
    ref = ref.upper()
    pool = [b for b in "ACGT" if b != ref]
    return rng.choice(pool)

def plot_obs_vs_sim(sim_counts, obs_counts, out_png):
    if plt is None:
        sys.stderr.write("[plot] matplotlib not available; skip plots.\n")
        return
    x = list(sim_counts)
    y = list(obs_counts)
    r_txt = "NA"; p_txt = "NA"
    if len(x) >= 2 and len(x) == len(y):
        if pearsonr is not None:
            r, p = pearsonr(x, y)
            r_txt = f"{r:.3f}"; p_txt = f"{p:.2e}"
        else:
            if np is not None:
                r = np.corrcoef(np.asarray(x), np.asarray(y))[0,1]
            else:
                n = len(x)
                mx = sum(x)/n; my = sum(y)/n
                num = sum((xi-mx)*(yi-my) for xi,yi in zip(x,y))
                denx = math.sqrt(sum((xi-mx)**2 for xi in x))
                deny = math.sqrt(sum((yi-my)**2 for yi in y))
                r = num/(denx*deny) if denx>0 and deny>0 else float('nan')
            r_txt = f"{r:.3f}"; p_txt = "NA"
            sys.stderr.write("[plot] scipy not available; Pearson p-value set to NA.\n")
    plt.figure()
    plt.scatter(x, y, s=8)
    plt.xlabel("Simulated variants per window")
    plt.ylabel("Observed variants per window")
    plt.title(f"Observed vs Simulated per-window counts\nPearson r={r_txt}, p={p_txt}")
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close()

def plot_sim_density_by_chrom(per_chrom_windows, sim_counts_by_chrom, out_prefix):
    if plt is None:
        return
    for chrom, wlist in per_chrom_windows.items():
        sim_counts = sim_counts_by_chrom.get(chrom, [0]*len(wlist))
        xs = []; ys = []
        for i, w in enumerate(wlist):
            call = max(1, w.callable)
            dens = float(sim_counts[i]) / call
            mid = (w.start + w.end) / 2.0
            xs.append(mid); ys.append(dens)
        plt.figure()
        plt.scatter(xs, ys, s=6)
        plt.xlabel("Position on chromosome (bp)")
        plt.ylabel("Simulated variant density (per bp)")
        plt.title(f"Simulated density by window: {chrom}")
        plt.tight_layout()
        out_png = f"{out_prefix}.sim_density.{chrom}.png"
        plt.savefig(out_png, dpi=160)
        plt.close()

def main():
    ap = argparse.ArgumentParser(
        description="Simulate heterozygous SNPs matching real per-window heterogeneity via inhomogeneous Poisson; support no-variant regions, optional mutated subset, & QC plots."
    )
    ap.add_argument("--ref", required=True, help="Reference FASTA (.fa/.fa.gz)")
    ap.add_argument("--win", type=int, default=100000, help="Window size (bp)")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--obs-vcf", help="Real SNPs VCF(.gz) to learn heterogeneity")
    src.add_argument("--obs-tsv", help="TSV: chrom start end count (per window)")
    tgt = ap.add_mutually_exclusive_group(required=True)
    tgt.add_argument("--heterozygosity", type=float, help="Target heterozygosity per bp (e.g. 0.01)")
    tgt.add_argument("--total-snps", type=int, help="Target total SNP count R")
    ap.add_argument("--alpha", type=float, default=0.0, help="Pseudocount smoothing for window weights")
    ap.add_argument("--exact-total", action="store_true", help="Use conditional Multinomial to enforce exact total R")
    ap.add_argument("--poisson", action="store_true", help="Use independent Poisson sampling (total fluctuates)")
    ap.add_argument("--seed", type=int, default=None, help="Random seed (default: None)")
    ap.add_argument("--hap1", default="hap1.fa")
    ap.add_argument("--hap2", default="hap2.fa")
    ap.add_argument("--vcf",  default="variants.vcf")  # .gz 可自动写 gz
    ap.add_argument("--plot-prefix", default=None, help="If set, write QC plots with this prefix")
    ap.add_argument("--nonvar-total", type=int, default=0, help="Total number of random no-variant regions to add (default 0)")
    ap.add_argument("--nonvar-bed",   default="nonvariant_regions.bed", help="Output BED for all no-variant regions")
    mutgrp = ap.add_mutually_exclusive_group()
    mutgrp.add_argument("--nonvar-mutate-frac", type=float, default=0.0, help="Fraction (0-1) of no-variant regions to receive 1–2 SNPs")
    mutgrp.add_argument("--nonvar-mutate-count", type=int, default=0, help="Absolute number of no-variant regions to receive 1–2 SNPs")
    ap.add_argument("--mut-nonvar-bed", default="mutated_nonvariant_regions.bed", help="BED of mutated no-variant regions (chrom start end length mut_count)")
    ap.add_argument("--mut-nonvar-variants-bed", default="mutated_nonvariant_variants.bed", help="BED of SNPs placed inside mutated no-variant regions (chrom start end ref>alt)")
    ap.add_argument("--out-count-dist", default="win_count_distribution.tsv",
                    help="Output TSV: per-window simulated variant count distribution (count<TAB>num_windows)")
    args = ap.parse_args()

    if args.poisson and args.exact_total:
        sys.stderr.write("Choose either --poisson or --exact-total, not both.\n"); sys.exit(2)
    if (not args.poisson) and (not args.exact_total):
        args.exact_total = True

    rng = random.Random(args.seed)
    if np is not None:
        np.random.seed(args.seed)

    order, seqs = read_fasta(args.ref)
    per_chrom_windows, all_windows = make_windows(seqs, order, args.win)
    callable_total = sum(w.callable for w in all_windows)
    if callable_total == 0:
        sys.stderr.write("No callable A/C/G/T bases in reference.\n"); sys.exit(2)

    if args.obs_vcf:
        obs = load_observed_from_vcf(args.obs_vcf, per_chrom_windows, args.win)
    else:
        obs = load_observed_from_tsv(args.obs_tsv, per_chrom_windows, args.win)

    if args.total_snps is not None:
        R = int(args.total_snps)
    else:
        R = int(round(args.heterozygosity * callable_total))

    O = []; C = []; WIndex = []
    for chrom in order:
        lst = per_chrom_windows[chrom]
        oi = obs.get(chrom, [0]*len(lst))
        if len(oi) != len(lst):
            if len(oi) < len(lst):
                oi = oi + [0]*(len(lst)-len(oi))
            else:
                oi = oi[:len(lst)]
        for k, w in enumerate(lst):
            O.append(max(0.0, float(oi[k])))
            C.append(float(w.callable))
            WIndex.append((chrom, k))
    sumO = sum(O)
    if sumO == 0:
        sys.stderr.write("Observed counts are all zero. Provide a real VCF/TSV or adjust inputs.\n"); sys.exit(2)

    if args.alpha > 0:
        meanC = (sum(C)/len(C)) if C else 1.0
        W = [O[i] + args.alpha * (C[i]/meanC) for i in range(len(O))]
    else:
        W = O[:]
    sumW = sum(W)
    if sumW <= 0:
        sys.stderr.write("All window weights are zero after smoothing.\n"); sys.exit(2)

    length_choices = [500] + list(range(1500, 20001, 1000))  # 500..20000 step 1000
    nonvar_regions = place_nonvariant_regions(
        seqs, order, args.nonvar_total, rng, length_choices
    )
    mask_index = build_interval_index(nonvar_regions)

    lambdas = [R * (w/sumW) for w in W]
    if args.poisson:
        K = poisson_allocate(lambdas, rng)
        R_alloc = sum(K)
    else:
        probs = [lam/sum(lambdas) for lam in lambdas]
        K = multinomial_allocate(R, probs)
        R_alloc = sum(K)

    hap1 = {chrom: bytearray(seqs[chrom].upper().encode("ascii")) for chrom in order}
    hap2 = {chrom: bytearray(seqs[chrom].upper().encode("ascii")) for chrom in order}
    variants = []  # (chrom, pos0, ref, alt)

    global_idx_of = { (chrom, k): idx for idx, (chrom, k) in enumerate(WIndex) }
    sim_counts_global = [0]*len(WIndex)
    sim_counts_by_chrom = {chrom: [0]*len(per_chrom_windows[chrom]) for chrom in order}

    for idx, k in enumerate(K):
        if k <= 0:
            continue
        chrom, k_local = WIndex[idx]
        w = per_chrom_windows[chrom][k_local]
        pos_list = choose_positions_excluding(
            seqs[chrom], w.start, w.end, k, rng, mask_index[chrom]
        )
        for pos0 in pos_list:
            ref = chr(hap1[chrom][pos0])
            if not is_callable_base(ref):
                continue
            alt = random_alt_base(ref, rng)
            if rng.random() < 0.5:
                hap1[chrom][pos0] = ord(alt)
            else:
                hap2[chrom][pos0] = ord(alt)
            variants.append((chrom, pos0, ref.upper(), alt))
            gi = global_idx_of[(chrom, k_local)]
            sim_counts_global[gi] += 1
            sim_counts_by_chrom[chrom][k_local] += 1

    all_nonvar_flat = []
    for chrom in order:
        for s,e,L in nonvar_regions.get(chrom, []):
            all_nonvar_flat.append( (chrom, s, e, L) )
    total_nonvar = len(all_nonvar_flat)
    want_mut = 0
    if total_nonvar > 0:
        if args.nonvar_mutate_count > 0:
            want_mut = min(args.nonvar_mutate_count, total_nonvar)
        elif args.nonvar_mutate_frac > 0.0:
            want_mut = min(total_nonvar, int(round(args.nonvar_mutate_frac * total_nonvar)))
    mutated_subset = []
    if want_mut > 0:
        sel = rng.sample(range(total_nonvar), want_mut)
        for idx in sel:
            mutated_subset.append(all_nonvar_flat[idx])

    mutated_regions_records = []  # lines for BED (chrom, s, e, L, mut_count)
    mutated_variant_records = []  # lines for BED (chrom, start, end, ref>alt)

    for chrom, s, e, L in mutated_subset:
        # 在该区间内放 1 或 2 个
        m = 1 if rng.random() < 0.5 else 2
        pos_list = choose_positions_in_interval(seqs[chrom], s, e, m, rng)
        mut_done = 0
        for pos0 in pos_list:
            ref = chr(hap1[chrom][pos0])
            if not is_callable_base(ref):
                continue
            alt = random_alt_base(ref, rng)
            if rng.random() < 0.5:
                hap1[chrom][pos0] = ord(alt)
            else:
                hap2[chrom][pos0] = ord(alt)
            variants.append((chrom, pos0, ref.upper(), alt))
            k_local = min(pos0 // args.win, len(per_chrom_windows[chrom]) - 1)
            gi = global_idx_of[(chrom, k_local)]
            sim_counts_global[gi] += 1
            sim_counts_by_chrom[chrom][k_local] += 1
            mutated_variant_records.append( (chrom, pos0, pos0+1, f"{ref.upper()}>{alt}") )
            mut_done += 1
        mutated_regions_records.append( (chrom, s, e, L, mut_done) )

    variants.sort(key=lambda x: (x[0], x[1]))

    def write_fa(path, seqdict):
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "wt") as out:
            for chrom in order:
                out.write(f">{chrom}\n")
                s = seqdict[chrom].decode("ascii")
                for i in range(0, len(s), 60):
                    out.write(s[i:i+60] + "\n")

    write_fa(args.hap1, hap1)
    write_fa(args.hap2, hap2)

    is_gz = args.vcf.endswith(".gz")
    opener = gzip.open if is_gz else open
    with opener(args.vcf, "wt") as v:
        v.write("##fileformat=VCFv4.2\n")
        v.write("##source=inhomogeneous_poisson_window_sim\n")
        v.write("##INFO=<ID=WIN,Number=3,Type=Integer,Description=\"window_start,window_end,assigned_K\">\n")
        v.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        for chrom in order:
            v.write(f"##contig=<ID={chrom},length={len(seqs[chrom])}>\n")
        v.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIM\n")

        K_by_window = defaultdict(int)
        for idx, k in enumerate(K):
            if k>0:
                chrom, k_local = WIndex[idx]
                K_by_window[(chrom, k_local)] = k

        for chrom, pos0, ref, alt in variants:
            k_local = min(pos0 // args.win, len(per_chrom_windows[chrom]) - 1)
            win = per_chrom_windows[chrom][k_local]
            Ki = K_by_window.get((chrom, k_local), 0)
            pos1 = pos0 + 1
            v.write(f"{chrom}\t{pos1}\t.\t{ref}\t{alt}\t.\tPASS\tWIN={win.start},{win.end},{Ki}\tGT\t0/1\n")

    if args.nonvar_total > 0:
        with open(args.nonvar_bed, "w") as bed:
            for chrom in order:
                for s,e,L in nonvar_regions.get(chrom, []):
                    bed.write(f"{chrom}\t{s}\t{e}\t{L}\n")

    if mutated_regions_records:
        with open(args.mut_nonvar_bed, "w") as bed:
            for chrom, s, e, L, mcount in mutated_regions_records:
                bed.write(f"{chrom}\t{s}\t{e}\t{L}\t{mcount}\n")
    if mutated_variant_records:
        with open(args.mut_nonvar_variants_bed, "w") as bedv:
            for chrom, s, e, tag in mutated_variant_records:
                bedv.write(f"{chrom}\t{s}\t{e}\t{tag}\n")

    if args.plot_prefix is not None:
        sim_counts = sim_counts_global
        obs_counts = [float(x) for x in O]
        out1 = f"{args.plot_prefix}.obs_vs_sim.png"
        plot_obs_vs_sim(sim_counts, obs_counts, out1)
        plot_sim_density_by_chrom(per_chrom_windows, sim_counts_by_chrom, args.plot_prefix)

    dist = Counter(sim_counts_global)
    with open(args.out_count_dist, "w") as fdist:
        for k in sorted(dist.keys()):
            fdist.write(f"{k}\t{dist[k]}\n")

    simulated = len(variants)
    sys.stderr.write(
        f"[done] callable={callable_total}  R_target={R}  allocated={R_alloc}  simulated={simulated}  "
        f"no_variant_regions={args.nonvar_total}  mutated_regions={len(mutated_regions_records)}  "
        f"windows={len(all_windows)}  dist_out={args.out_count_dist}\n"
    )

if __name__ == "__main__":
    main()
