#!/bin/bash
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=pipeline2
#SBATCH --output=log/pipeline2_%j.out


prefix=Chr${SLURM_ARRAY_TASK_ID}
[ "$prefix" = "Chr16" ] && prefix="Chr1620"

# 并行执行染色体处理
python 3_Concatenation_CDS_all_individuals.py --part ${SLURM_ARRAY_TASK_ID}

# cat results/CDS/CDS*_reverse.fasta > results/CDS/all_CDS_reverse.fasta
# cat results/CDS/CDS*_forward.fasta > results/CDS/all_CDS_forward.fasta