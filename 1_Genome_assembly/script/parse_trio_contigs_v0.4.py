#!/usr/bin/env python3
import pysam # type: ignore
import sys
import os
import tempfile
import argparse
from multiprocessing import Pool
from collections import defaultdict

# Version 0.4 changes:
# 1. Added processing of hard clipping in find_contig_position function

def parse_arguments():
    """Parse command line arguments for the script"""
    parser = argparse.ArgumentParser(
        description='Process trio contigs with parallel processing support',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--input', '-i', required=True,
                      help='Input BAM file path')
    parser.add_argument('--output', '-o', required=True,
                      help='Main output text file path')
    parser.add_argument('--uncallable', '-u', required=True,
                      help='Uncallable positions output path')
    parser.add_argument('--reference', '-r', required=True,
                      help='Reference genome FASTA file path')
    parser.add_argument('--threads', '-t', type=int, default=os.cpu_count(),
                      help='Number of parallel processing threads')
    parser.add_argument('--tmp-dir', default=tempfile.gettempdir(),
                      help='Temporary directory path')
    return parser.parse_args()

def get_parent_from_contig(contig_id: str):
    """Determine parent origin based on contig ID suffix"""
    suffix = contig_id[-4:]
    return "mother" if suffix == "Bf_M" else "father" if suffix == "Bf_P" else None

def find_contig_position(read, target_ref_pos):
    """Find corresponding position in contig for given reference position"""
    if read.is_unmapped:
        return None, None

    ref_pos = read.reference_start
    contig_pos = 0
    seq_pos = 0

    # Process CIGAR operations to map reference position to contig position
    for op, length in read.cigartuples:
        if ref_pos > target_ref_pos:
            break
            
        if op in (0, 7, 8):  # Match/mismatch operations
            op_end = ref_pos + length
            if ref_pos <= target_ref_pos < op_end:
                offset = target_ref_pos - ref_pos
                final_contig_pos = contig_pos + offset
                final_seq_pos = seq_pos + offset
                return final_contig_pos, read.query_sequence[final_seq_pos]
            ref_pos += length
            contig_pos += length
            seq_pos += length
        elif op in (1, 4):  # Insertion or soft-clip
            contig_pos += length
            seq_pos += length
        elif op in (2, 3):  # Deletion or reference skip
            ref_pos += length
        elif op == 5:  # Hard clip
            contig_pos += length  # Only affects contig coordinates
    return None, None

def process_region(args):
    """Process a single genomic region in an independent worker process"""
    bam_path, ref_path, chrom, start, end, output_tmp, uncallable_tmp = args
    try:
        with (pysam.AlignmentFile(bam_path, "rb") as samfile,
             pysam.FastaFile(ref_path) as ref_fasta,
             open(output_tmp, 'w') as out_fd,
             open(uncallable_tmp, 'w') as uncall_fd):

            # Iterate through all positions in the specified region
            for pileup_col in samfile.pileup(chrom, start, end, stepper="nofilter"):
                ref_pos_0based = pileup_col.reference_pos
                if not (start <= ref_pos_0based < end):
                    continue

                ref_pos_1based = ref_pos_0based + 1
                try:
                    ref_base = ref_fasta.fetch(chrom, ref_pos_0based, ref_pos_0based+1).upper()
                except KeyError:
                    ref_base = "N"

                parent_map = defaultdict(list)
                for pileup_read in pileup_col.pileups:
                    if pileup_read.is_del or pileup_read.is_refskip:
                        continue
                    
                    read = pileup_read.alignment
                    parent = get_parent_from_contig(read.query_name)
                    if not parent:
                        continue
                    
                    contig_pos, base = find_contig_position(read, ref_pos_0based)
                    if contig_pos is None or base is None:
                        continue
                    
                    # Handle reverse strand orientation
                    if read.is_reverse:
                        contig_length = read.infer_query_length(always=True) 
                        contig_pos = contig_length - contig_pos - 1
                    
                    parent_map[parent].append((
                        read.query_name,
                        contig_pos + 1,  # Convert to 1-based position
                        base
                    ))

                # Determine if position is callable
                mother_count = len(parent_map.get("mother", []))
                father_count = len(parent_map.get("father", []))
                total = mother_count + father_count
                uncallable = False
                
                # Uncallable conditions based on read counts
                if total > 4 or (total == 3 and (mother_count == 3 or father_count == 3)) \
                    or (total == 4 and (mother_count > 2 or father_count > 2)) \
                    or (total == 2 and (mother_count > 1 or father_count > 1)) or total == 0:
                    uncall_fd.write(f"{chrom}\t{ref_pos_1based}\n")
                    continue
                
                # Format output line with allele information
                def fmt_alleles(alleles):
                    return [f"{cid}:{pos}:{base}" for cid, pos, base in alleles[:2]] + ["NA"]*(2-len(alleles))
                
                out_line = "\t".join([
                    chrom,
                    str(ref_pos_1based),
                    ref_base,
                    *fmt_alleles(parent_map.get("mother", [])),
                    *fmt_alleles(parent_map.get("father", [])),
                    ref_base
                ]) + "\n"
                out_fd.write(out_line)

        return True
    except Exception as e:
        print(f"Error processing {chrom}:{start}-{end}: {str(e)}")
        return False

def main():
    args = parse_arguments()
    
    # Ensure temporary directory exists
    os.makedirs(args.tmp_dir, exist_ok=True)

    # Get chromosome length information from BAM header
    with pysam.AlignmentFile(args.input, "rb") as samfile:
        chrom_info = [(name, length) for name, length in zip(samfile.references, samfile.lengths)]
    
    # Generate parallel task parameters (1MB chunks)
    TASK_CHUNK = 1_000_000
    tasks = []
    tmp_files = []
    
    # Create temporary files in specified temp directory
    for chrom, length in chrom_info:
        for chunk_start in range(0, length, TASK_CHUNK):
            chunk_end = min(chunk_start + TASK_CHUNK, length)
            output_tmp = tempfile.mktemp(
                suffix=f".{chrom}_{chunk_start}.out", 
                dir=args.tmp_dir
            )
            uncallable_tmp = tempfile.mktemp(
                suffix=f".{chrom}_{chunk_start}.uncall",
                dir=args.tmp_dir
            )
            tasks.append((args.input, args.reference, chrom, chunk_start, chunk_end, output_tmp, uncallable_tmp))
            tmp_files.append((output_tmp, uncallable_tmp))
    
    # Parallel processing with worker pool
    with Pool(processes=args.threads) as pool:
        results = pool.map(process_region, tasks)
    
    # Merge results from temporary files
    with open(args.output, 'w') as out_final, open(args.uncallable, 'w') as uncall_final:
        # Write output header
        out_final.write("\t".join([
            "ref_chr", "ref_pos", "ref_allele",
            "mother_allele1", "mother_allele2",
            "father_allele1", "father_allele2",
            "genotype"]) + "\n")
        
        # Merge files in chromosomal order
        for output_tmp, uncallable_tmp in tmp_files:
            if os.path.exists(output_tmp):
                with open(output_tmp) as f:
                    out_final.write(f.read())
                os.remove(output_tmp)
            if os.path.exists(uncallable_tmp):
                with open(uncallable_tmp) as f:
                    uncall_final.write(f.read())
                os.remove(uncallable_tmp)

if __name__ == "__main__":
    main()
