#!/bin/bash
#SBATCH -w node80
#SBATCH -c 2
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=freebayes
#SBATCH --output=log/freebayes_%j_%a.out


python callable_genome.py Bf-${SLURM_ARRAY_TASK_ID}    
bedtools intersect \
    -a /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/FNR/callable_genome/bed/Bf-${SLURM_ARRAY_TASK_ID}.contig_inheritance_v0.8.postprocess.remove_similarity.bed \
    -b /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/default/callable_genome/progeny_depth/result/Bf-${SLURM_ARRAY_TASK_ID}.contig_inheritance_v0.8.postprocess.remove3foldsites.sort.bed |\
    sort -k1,1 -k2,2n > bed/Bf-${SLURM_ARRAY_TASK_ID}.contig_inheritance_v0.8.postprocess.remove_similarity.remove3foldsites.bed