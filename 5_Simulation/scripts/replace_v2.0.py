# -*- coding: utf-8 -*-
import random
import argparse
import sys
import gzip
import os
from collections import defaultdict

def read_fasta(file_path):
    """Read FASTA file and return a dictionary of sequence IDs and sequences."""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]  # Take only the first word as ID
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
    return sequences

def read_bed(bed_file):
    """
    Read BED file and return dictionary of positions per sequence ID.
    Format: {seq_id: {position1: True, position2: True, ...}}
    """
    bed_positions = defaultdict(set)
    
    opener = gzip.open if bed_file.endswith('.gz') else open
    with opener(bed_file, 'rt') as f:  # 'rt' mode for text reading
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue  # Skip malformed lines
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            
            # Handle different BED formats (0-based, half-open)
            for pos in range(start, end):
                bed_positions[chrom].add(pos)
                
    return bed_positions

def calculate_genome_mutations(sequences, bed_positions, total_mutations):
    all_valid_positions = []
    for seq_id, seq in sequences.items():
        valid_positions = bed_positions.get(seq_id, set())
        available_positions = [pos for pos in valid_positions if pos < len(seq)]
        all_valid_positions.extend([(seq_id, pos) for pos in available_positions])
    
    total_available = len(all_valid_positions)
    if total_available == 0:
        raise ValueError("No valid mutation positions found in BED regions")
    if total_mutations > total_available:
        raise ValueError(f"Requested {total_mutations} mutations but only {total_available} positions available in BED regions across all sequences")
    
    selected_positions = random.sample(all_valid_positions, total_mutations)
    
    mutations_per_seq = defaultdict(list)
    for seq_id, pos in selected_positions:
        mutations_per_seq[seq_id].append(pos)
    
    return mutations_per_seq

def apply_mutations(sequences, mutation_assignments):
    mutated_sequences = sequences.copy()
    mutation_records = defaultdict(list)
    bases = ['A', 'T', 'C', 'G']
    
    for seq_id, positions in mutation_assignments.items():
        if seq_id not in sequences:
            continue
            
        seq_list = list(sequences[seq_id])  
        
        for pos in positions:
            if pos >= len(seq_list):
                continue
                
            original_char = seq_list[pos]
            original_base = original_char.upper()
            
            choices = [b for b in bases if b != original_base]
            
            if choices:
                new_base = random.choice(choices)
                seq_list[pos] = new_base if original_char.isupper() else new_base.lower()
                mutation_records[seq_id].append((pos, original_base, new_base))
        
        mutated_sequences[seq_id] = ''.join(seq_list)
    
    return mutated_sequences, mutation_records

def write_fasta(sequences, output_file):
    """Write sequences to FASTA file with 80-character line width."""
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    with open(output_file, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f'>{seq_id}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

def write_mutations(mutation_records, output_file):
    """Write mutation records to a tab-separated file."""
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    with open(output_file, 'w') as f:
        f.write("Sequence\tPosition\tOriginal_Base\tNew_Base\n")
        for seq_id, mutations in mutation_records.items():
            for pos, original, new in mutations:
                f.write(f"{seq_id}\t{pos+1}\t{original}\t{new}\n")

def main():
    parser = argparse.ArgumentParser(description='Introduce specified number of base substitutions in BED regions')
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('bed_file', help='BED file defining mutation regions')
    parser.add_argument('num_mutations', type=int, help='Total number of mutations to introduce in the genome')
    parser.add_argument('-p', '--prefix', default='mutated_sequence', 
                        help='Output file prefix (default: mutated_sequence), will generate <prefix>.fa and <prefix>.mut.txt')
    
    args = parser.parse_args()
    
    try:
        sequences = read_fasta(args.input_file)
        if not sequences:
            sys.exit("Error: No valid sequences found in input file")
            
        bed_positions = read_bed(args.bed_file)
        
        mutation_assignments = calculate_genome_mutations(
            sequences, 
            bed_positions, 
            args.num_mutations
        )
        
        mutated_sequences, mutation_records = apply_mutations(
            sequences,
            mutation_assignments
        )
        
        fasta_output = f"{args.prefix}.fa"
        report_output = f"{args.prefix}.mut.txt"
        
        write_fasta(mutated_sequences, fasta_output)
        write_mutations(mutation_records, report_output)
        
        actual_mutations = sum(len(muts) for muts in mutation_records.values())
        mutated_chromosomes = len(mutation_records)
        total_chromosomes = len(sequences)
        
        print(f"Successfully introduced {actual_mutations} mutations across {mutated_chromosomes}/{total_chromosomes} sequences")
        print(f"Mutated FASTA: {fasta_output}")
        print(f"Mutation report: {report_output}")
    
    except Exception as e:
        sys.exit(f"Error: {str(e)}")

if __name__ == '__main__':
    main()
