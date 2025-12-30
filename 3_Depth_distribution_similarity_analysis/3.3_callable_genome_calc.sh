#!/bin/bash
#SBATCH -w node80
#SBATCH -c 2
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=freebayes
#SBATCH --output=log/freebayes_%j_%a.out


python scripts/callable_genome.v0.3.py \
  -i postprocess_result/${prefix}/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.out \
  -o postprocess_result/${prefix}/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3 \
  --parental-homology test/mat_pat_rbh.contig_pairs.list --len-col length 

bedtools subtract -a postprocess_result/${prefix}/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.callable_filtered.bed \
    -b /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/default/callable_genome/progeny_depth/tmp/${prefix}.tmp.bed > postprocess_result/${prefix}/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.callable_filtered.remove3foldsites.bed
