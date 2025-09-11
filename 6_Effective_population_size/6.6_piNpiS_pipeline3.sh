#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=pipeline3
#SBATCH --output=log/pipeline3_%j.out

prefix=Chr${SLURM_ARRAY_TASK_ID}

if [ "$prefix" = "Chr16" ]; then
    prefix="Chr1620"
fi
indi=$(head -n ${SLURM_ARRAY_TASK_ID} indi.list | tail -1)

python 4_catcon_chr.py --chr ${prefix}

cat results/CDS/${prefix}_forward.fasta results/CDS/${prefix}_reverse.fasta > results/CDS/${prefix}_forward_reverse.fasta

python 5_remove_N.py --chr ${prefix}

./dNdSpiNpiS_1.0 -alignment_file=results/CDS/${prefix}_forward_reverse_no_indiv_with_N.fasta -ingroup=bf results/${prefix}_forward_reverse_no_indiv_with_N.out results/${prefix}_forward_reverse_no_indiv_with_N.sum

python 6_Get_corrector_factor.py
