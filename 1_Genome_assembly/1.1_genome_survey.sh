#!/bin/bash
#SBATCH -c 20
#SBATCH --time=48:00:00                     
#SBATCH --mem=20G
#SBATCH --output=log/kmer_%j.out

prefix=$1
mkdir tmp
fastq=$2

# Count k-mers using Jellyfish (21-mer size)
jellyfish count -C -m 21 -s 50G -t 20 -o ${prefix}_reads.jf <(zcat fastq/${prefix}_qc2_1.fq.gz) <(zcat fastq/${prefix}_qc2_2.fq.gz)
# Generate k-mer histogram from Jellyfish results
jellyfish histo -t 20 ${prefix}_reads.jf > ${prefix}_reads.histo
# Run GenomeScope2 for genome analysis
genomescope2 -i ${prefix}_reads.histo -o ${prefix}_output_dir -k 21 