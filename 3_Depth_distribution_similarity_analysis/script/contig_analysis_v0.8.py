import re
import os
import gzip
import hashlib
import numpy as np # type: ignore
from collections import defaultdict
import bisect
import tempfile
from multiprocessing import Pool
import json
import argparse
from scipy import stats  # type: ignore
#1. add zero-depth contig.

class FilterLogger:
    def __init__(self):
        self.filtered_contigs = []
        self.total_filtered_length = 0
    
    def log_filtered(self, contig, length):
        self.filtered_contigs.append(contig)
        self.total_filtered_length += length
    
    def write_log(self, log_file):
        with open(log_file, 'w') as f:
            f.write(f"Filtered contigs: {len(self.filtered_contigs)}\n")
            f.write(f"Total filtered length: {self.total_filtered_length}\n")
            f.write("\nFiltered contig list:\n")
            f.write("\n".join(self.filtered_contigs))

def get_all_contigs(*temp_dirs):
    all_contigs = set()
    for temp_dir in temp_dirs:
        try:
            with open(os.path.join(temp_dir, 'metadata.json')) as f:
                all_contigs.update(json.load(f).values())
        except FileNotFoundError:
            continue
    return list(all_contigs)

def safe_contig_name(contig):
    """Generate safe filename while preserving original contig mapping"""
    return hashlib.md5(contig.encode()).hexdigest()

def split_bed_to_temp(input_file, temp_dir):
    """Split BED file into contig-based temporary files and save metadata"""
    os.makedirs(temp_dir, exist_ok=True)
    metadata = {}
    open_func = gzip.open if input_file.endswith('.gz') else open
    mode = 'rt' if input_file.endswith('.gz') else 'r'
    
    with open_func(input_file, mode) as f_in:
        current_contig = None
        current_file = None
        for line in f_in:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
            contig = parts[0]
            if contig != current_contig:
                if current_file:
                    current_file.close()
                current_contig = contig
                safe_name = safe_contig_name(contig)
                metadata[safe_name] = contig
                file_path = os.path.join(temp_dir, safe_name)
                current_file = open(file_path, 'a')
                current_file.write(f"#CONTIG={contig}\n")  
            current_file.write(line)
        if current_file:
            current_file.close()
    
    with open(os.path.join(temp_dir, 'metadata.json'), 'w') as f:
        json.dump(metadata, f)
    return temp_dir

def read_contig_data(temp_dir, safe_name):
    file_path = os.path.join(temp_dir, safe_name)
    contig = None
    data = []
    try:
        with open(file_path, 'r') as f:
            header = f.readline().strip()
            if header.startswith("#CONTIG="):
                contig = header.split('=')[1]
            for line in f:
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                s, e, depth = int(parts[1]), int(parts[2]), float(parts[3])
                data.append((s, e, depth))
    except FileNotFoundError:
        pass
    return contig, data

def get_common_contigs(P_temp, M_temp, pr_temp):
    def load_metadata(temp_dir):
        with open(os.path.join(temp_dir, 'metadata.json'), 'r') as f:
            return json.load(f)
    
    p_meta = load_metadata(P_temp)
    m_meta = load_metadata(M_temp)
    pr_meta = load_metadata(pr_temp)
    
    common = set(p_meta.keys()) & set(m_meta.keys()) & set(pr_meta.keys())
    return [(p_meta[s], s) for s in common if 
            p_meta[s] == m_meta[s] == pr_meta[s]]

def get_progeny_avg_depth(pr_id):
    file_path = f"/public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/default/bam/{pr_id}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam.depth.mosdepth.summary.txt"
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith("total"):
                    return float(line.split("\t")[3])
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    return None

def process_contig(args):
    orig_name, safe_name, pr_id, P_temp, M_temp, pr_temp, threshold, method  = args

    contig_length = get_contig_length(orig_name)
    if contig_length is None or contig_length < 500:
        return (orig_name, None)  
    
    def safe_read_contig_data(temp_dir, safe_name):
        try:
            return read_contig_data(temp_dir, safe_name)
        except FileNotFoundError:
            return None, []
    
    _, P_data = read_contig_data(P_temp, safe_name)
    _, M_data = read_contig_data(M_temp, safe_name)
    _, pr_data = read_contig_data(pr_temp, safe_name)
    
    _, P_data = safe_read_contig_data(P_temp, safe_name) or (None, [])
    _, M_data = safe_read_contig_data(M_temp, safe_name) or (None, [])
    _, pr_data = safe_read_contig_data(pr_temp, safe_name) or (None, [])
    
    def _process(contig, intervals):
        windows = split_contig_into_windows(get_contig_length(contig))
        if not windows:
            return None, 0
        if not intervals:
            return [0.0] * len(windows), 0.0
        win_depths = [calculate_interval_depth(intervals, s, e) for s, e in windows]
        total = sum(d*(e-s) for s,e,d in intervals)
        avg = total / contig_length if contig_length else 0
        return [d/avg for d in win_depths], avg
    
    P_norm, P_avg = _process(orig_name, P_data)
    M_norm, M_avg = _process(orig_name, M_data)
    pr_norm, pr_avg = _process(orig_name, pr_data)

    metrics = calculate_metrics(orig_name, pr_id, method, P_norm, M_norm, pr_norm, pr_avg, threshold)
    return (orig_name, metrics)

def calculate_metrics(contig, pr_id, method, P_norm, M_norm, pr_norm, pr_avg, threshold):
    progeny_avg_depth = get_progeny_avg_depth(pr_id)

    metrics = {}
    metrics['low_depth_num'] = 0
    metrics['low_depth_length'] = 0
    metrics['normal_depth_num'] = 0
    metrics['normal_depth_length'] = 0
    metrics['pass_num'] = 0
    metrics['pass_length'] = 0
    metrics['mother_pass_num'] = 0
    metrics['mother_pass_length'] = 0
    metrics['father_pass_num'] = 0
    metrics['father_pass_length'] = 0

    if float(pr_avg) < (progeny_avg_depth / 4) or float(pr_avg) < (progeny_avg_depth / 4) or float(pr_avg) < 10:
        depth = "low_depth"
        metrics['low_depth_num'] = 1
        metrics['low_depth_length'] = get_contig_length(contig)
    else:
        depth = "normal_depth"
        metrics['normal_depth_num'] = 1
        metrics['normal_depth_length'] = get_contig_length(contig)
    
    metrics['P_pearson'], metrics['P_pearson_p'] = safe_pearson(P_norm, pr_norm)
    metrics['M_pearson'], metrics['M_pearson_p'] = safe_pearson(M_norm, pr_norm)
    metrics['M_cosine'] = safe_cosine(M_norm, pr_norm)
    metrics['P_cosine'] = safe_cosine(P_norm, pr_norm)
    metrics['P_spearman'], metrics['P_spearman_p'] = safe_spearman(P_norm, pr_norm)
    metrics['M_spearman'], metrics['M_spearman_p'] = safe_spearman(M_norm, pr_norm)
    metrics['P_hybrid'] = zero_aware_similarity(P_norm, pr_norm)
    metrics['M_hybrid'] = zero_aware_similarity(M_norm, pr_norm)
    metrics['MP_pearson'], metrics['MP_pearson_p'] = safe_pearson(M_norm, P_norm)
    metrics['MP_hybrid'] = zero_aware_similarity(M_norm, P_norm)
    metrics['depth'] = depth
    metrics['length'] = get_contig_length(contig)
    metrics['parenthood'] = get_contig_parenthood(contig)

    m_chosen_method = f'M_{method}'
    p_chosen_method = f'P_{method}'     
        
    if metrics[m_chosen_method] > metrics[p_chosen_method] and metrics['parenthood'] == 'mother':
        metrics['inheritance'] = "PASS"
        metrics['pass_num'] = 1
        metrics['mother_pass_num'] = 1
        metrics['pass_length'] = metrics['length']
        metrics['mother_pass_length'] = metrics['length']
    elif metrics[p_chosen_method] > metrics[m_chosen_method] and metrics['parenthood'] == 'father':
        metrics['inheritance'] = "PASS"
        metrics['pass_num'] = 1
        metrics['father_pass_num'] = 1
        metrics['pass_length'] = metrics['length']
        metrics['father_pass_length'] = metrics['length']
    elif metrics[p_chosen_method] > metrics[m_chosen_method] and metrics['parenthood'] == 'mother':
        metrics['inheritance'] = "mother_low"
    elif metrics[m_chosen_method] > metrics[p_chosen_method] and metrics['parenthood'] == 'father':
        metrics['inheritance'] = "father_low"
    else:
        metrics['inheritance'] = "other"

    if metrics['MP_hybrid'] > threshold:
        metrics['MP_similar'] = "parental_Similar"
    else:
        metrics['MP_similar'] = "parental_Dissimilar"

    if metrics['depth']=="low_depth":
        metrics['inheritance'] = metrics['depth']
        metrics['pass_num'] = 0  
        metrics['pass_length'] = 0
        metrics['mother_pass_num'] = 0
        metrics['mother_pass_length'] = 0
        metrics['father_pass_num'] = 0
        metrics['father_pass_length'] = 0
    
    return metrics

def get_contig_length(contig_name):
    match = re.search(r'len(\d+)', contig_name)
    return int(match.group(1)) if match else None

def get_contig_parenthood(contig_name):
    if contig_name.endswith('M'):
        return 'mother'
    elif contig_name.endswith('P'):
        return 'father'
    else:
        return None

def split_contig_into_windows(contig_length, window_size=50):
    if contig_length < window_size:
        return []
    num_windows = contig_length // window_size
    return [(i*window_size, (i+1)*window_size) for i in range(num_windows)]

def calculate_interval_depth(sorted_intervals, win_start, win_end):
    win_len = win_end - win_start
    if win_len <= 0:
        return 0.0
    
    # Extract pre-generated start/end lists
    starts = [s for s, e, d in sorted_intervals]
    ends = [e for s, e, d in sorted_intervals]
    depths = [d for s, e, d in sorted_intervals]
    
    total = 0.0
    
    # Determine check range using bisect
    # First possible overlapping interval: ends[i] > win_start
    left = bisect.bisect_right(ends, win_start)
    # Last possible overlapping interval: starts[i] < win_end
    right = bisect.bisect_left(starts, win_end)
    
    for i in range(left, right):
        s, e, d = starts[i], ends[i], depths[i]
        overlap_start = max(s, win_start)
        overlap_end = min(e, win_end)
        overlap_len = overlap_end - overlap_start
        
        if overlap_len > 0:
            total += d * overlap_len
    
    return total / win_len

def safe_pearson(a, b, return_p=True, nan_policy='omit'):
    a = np.asarray(a)
    b = np.asarray(b)
    
    if np.all(a == a[0]) or np.all(b == b[0]):
        if return_p:
            return 0.0, 1.0 
        return 0.0
    
    mask = np.isnan(a) | np.isnan(b)
    if nan_policy == 'omit':
        a = a[~mask]
        b = b[~mask]
    elif np.any(mask):
        if return_p:
            return 0.0, 1.0
        return 0.0
    
    n = len(a)
    if n < 2:
        if return_p:
            return 0.0, 1.0
        return 0.0
    
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.corrcoef(a, b)[0, 1]
    
    if return_p:
        if np.isnan(r) or abs(r) == 1.0:
            return 0.0, 1.0
        try:
            t_stat = r * np.sqrt((n - 2) / (1 - r**2))
            p = stats.t.sf(abs(t_stat), df=n-2)*2
            return np.nan_to_num(r, nan=0.0), np.nan_to_num(p, nan=1.0)
        except:
            return 0.0, 1.0
    return np.nan_to_num(r, nan=0.0)


def safe_cosine(a, b):
    norm_a, norm_b = np.linalg.norm(a), np.linalg.norm(b)
    if norm_a == 0 or norm_b == 0:
        return 0.0
    return np.dot(a, b) / (norm_a * norm_b)

# Newly added weighted Jaccard calculation function
def safe_weighted_jaccard(a, b):
    a = np.array(a)
    b = np.array(b)
    
    # Handle all-zero case
    if np.all(a == 0) and np.all(b == 0):
        return 1.0 
    
    min_sum = np.sum(np.minimum(a, b))
    max_sum = np.sum(np.maximum(a, b))
    
    if max_sum == 0:
        return 0.0
    return min_sum / max_sum

def zero_aware_similarity(a, b):
    a = np.array(a)
    b = np.array(b)
    
    zero_mask = (a == 0) | (b == 0)
    non_zero = ~zero_mask
    
    jaccard = safe_weighted_jaccard(a[zero_mask], b[zero_mask]) if zero_mask.any() else 1.0

    if non_zero.any():
        cosine = safe_cosine(a[non_zero], b[non_zero])
    else:
        cosine = 1.0
    w_j = np.sqrt(zero_mask.mean())
    w_c = 1 - w_j
    
    return w_j * jaccard + w_c * cosine

def safe_spearman(a, b):
    with np.errstate(invalid='ignore'):
        if len(a) < 3 or len(b) < 3:
            return 0.0, 1.0
        
        try:
            corr, p = stats.spearmanr(a, b, nan_policy='omit')
            corr = 0.0 if np.isnan(corr) else corr
            p = 1.0 if np.isnan(p) else p
            return corr, p
        except:
            return 0.0, 1.0


def write_results(output_file, stats_file,results):
    def format_value(value, fmt=".3f"):
        if value is None:
            return "NA"
        try:
            return f"{value:{fmt}}"
        except (TypeError, ValueError):
            return "NA"

    try:
        with open(output_file, 'w') as f:
            headers = [
                "ContigID",
                "P_Pearson", "P_Pearson_p",
                "M_Pearson", "M_Pearson_p",
                "P_Cosine","M_Cosine",
                "P_Spearman", "P_Spearman_p",
                "M_Spearman", "M_Spearman_p",
                "P_Hybrid", "M_Hybrid",
                "MP_Pearson", "MP_Pearson_p",
                "MP_hybrid",
                "length",
                "parenthood",
                "inheritance",
                "MP_Similar"
            ]
            f.write("\t".join(headers) + "\n")
            
            non_num = 0

            for contig, metrics in results:
                if metrics is None:
                    non_num+=1
                    continue

                fields = [
                    contig,
                    format_value(metrics.get('P_pearson')),
                    format_value(metrics.get('P_pearson_p')),
                    format_value(metrics.get('M_pearson')),
                    format_value(metrics.get('M_pearson_p')),
                    format_value(metrics.get('P_cosine')),
                    format_value(metrics.get('M_cosine')),
                    format_value(metrics.get('P_spearman')),
                    format_value(metrics.get('P_spearman_p')),
                    format_value(metrics.get('M_spearman')),
                    format_value(metrics.get('M_spearman_p')),
                    format_value(metrics.get('P_hybrid')),
                    format_value(metrics.get('M_hybrid')),
                    format_value(metrics.get('MP_pearson')),
                    format_value(metrics.get('MP_pearson_p')),
                    format_value(metrics.get('MP_hybrid')),
                    str(metrics.get('length', 'NA')),
                    str(metrics.get('parenthood', 'NA')),
                    str(metrics.get('inheritance', 'NA')),
                    str(metrics.get('MP_similar', 'NA'))
                ]
                f.write("\t".join(fields) + "\n")

        f.close()

        low_depth_num = 0
        low_depth_length = 0
        normal_depth_num = 0
        normal_depth_length = 0
        pass_num = 0
        pass_length = 0    
        mother_pass_num = 0
        mother_pass_length = 0
        father_pass_num = 0
        father_pass_length = 0

        for contig, metrics in results:
            if metrics is None:
                continue
            low_depth_num += metrics['low_depth_num']
            low_depth_length += metrics['low_depth_length']
            normal_depth_num += metrics['normal_depth_num']
            normal_depth_length += metrics['normal_depth_length']
            pass_num += metrics['pass_num'] 
            pass_length += metrics['pass_length']    
            mother_pass_num += metrics['mother_pass_num']
            mother_pass_length += metrics['mother_pass_length']
            father_pass_num += metrics['father_pass_num']
            father_pass_length += metrics['father_pass_length']

        with open(stats_file, 'w') as f:
            f.write(f"low_depth_num\t{low_depth_num}\n")
            f.write(f"low_depth_length\t{low_depth_length}\n")
            f.write(f"normal_depth_num\t{normal_depth_num}\n")
            f.write(f"normal_depth_length\t{normal_depth_length}\n")
            f.write(f"pass_num\t{pass_num}\n")
            f.write(f"pass_length\t{pass_length}\n")
            f.write(f"mother_pass_num\t{mother_pass_num}\n")
            f.write(f"mother_pass_length\t{mother_pass_length}\n")
            f.write(f"father_pass_num\t{father_pass_num}\n")
            f.write(f"father_pass_length\t{father_pass_length}\n")

        f.close()    
        return True
    
    except Exception as e:
        print(f"Unknown error occurred while writing file: {str(e)}")

    return False

def main(P, M, progeny, output, pr_id, method, threshold=0.9, processes=4):

    logger = FilterLogger()

    with tempfile.TemporaryDirectory(dir="/public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/cosin_similarity/tmp") as P_temp, \
         tempfile.TemporaryDirectory(dir="/public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/cosin_similarity/tmp") as M_temp, \
         tempfile.TemporaryDirectory(dir="/public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/cosin_similarity/tmp") as pr_temp:

        split_bed_to_temp(P, P_temp)
        split_bed_to_temp(M, M_temp)
        split_bed_to_temp(progeny, pr_temp)
        
        all_contigs = get_all_contigs(P_temp, M_temp, pr_temp)
        args = [(contig, safe_contig_name(contig), pr_id, P_temp, M_temp, pr_temp, threshold, method)
            for contig in all_contigs]
        

        with Pool(processes=processes) as pool:
            results = pool.map(process_contig, args)

        for orig_name, metrics in results:
            if metrics is None:
                contig_length = get_contig_length(orig_name) or 0
                logger.log_filtered(orig_name, contig_length)

        log_file = os.path.splitext(output)[0] + ".log"
        logger.write_log(log_file)
        stats_file = os.path.splitext(output)[0] + ".stats"
        
        if output:
            write_results(output, stats_file, results)
        else:
            print("ContigID\tP_Pearson\tM_Pearson\tP_Cosine\tM_Cosine\tP_WeightedJaccard\tM_WeightedJaccard\tP_Hybrid\tM_Hybrid\tMP_WeightedJaccard\tInheritance\tMP_Similair")
            for contig, metrics in results:
                print(f"{contig}\t"
                    f"{metrics.get('P_pearson', np.nan):.3f}",
                    f"{metrics.get('P_pearson_p', np.nan):.3f}",
                    f"{metrics.get('M_pearson', np.nan):.3f}",
                    f"{metrics.get('M_pearson_p', np.nan):.3f}", 
                    f"{metrics.get('P_spearman', np.nan):.3f}",
                    f"{metrics.get('P_spearman_p', np.nan):.3f}",
                    f"{metrics.get('M_spearman', np.nan):.3f}",
                    f"{metrics.get('M_spearman_p', np.nan):.3f}",
                    f"{metrics.get('P_hybrid', np.nan):.3f}",
                    f"{metrics.get('M_hybrid', np.nan):.3f}",
                    f"{metrics.get('MP_hybrid', np.nan):.3f}",
                    f"{metrics.get('MP_spearman', np.nan):.3f}",
                    f"{metrics.get('MP_spearman_p', np.nan):.3f}",
                    str(metrics.get('length', 'NA')),
                    str(metrics.get('parenthood', 'NA')),
                    str(metrics.get('inheritance', 'NA')),
                    str(metrics.get('MP_similair', 'NA')))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--parent", required=True)
    parser.add_argument("-m", "--mother", required=True)
    parser.add_argument("-pr", "--progeny", required=True)
    parser.add_argument("--pr_id", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--threshold", type=float, default=0.9)
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--method", choices=['pearson', 'spearman', 'hybrid'], default='hybrid')

    args = parser.parse_args()
    
    main(args.parent, args.mother, args.progeny, args.output, args.pr_id, args.method,
        args.threshold, processes=args.processes)
