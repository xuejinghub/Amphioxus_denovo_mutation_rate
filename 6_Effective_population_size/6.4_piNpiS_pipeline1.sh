#!/bin/bash
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=pipeline1
#SBATCH --output=log/pipeline1_%j.out

indi=$(head -n ${SLURM_ARRAY_TASK_ID} indi.list | tail -1)

# seqkit replace --ignore-case --kv-file fasta.rename.list --pattern "^(\w+)" --replacement "{kv}" /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/effective_population_size/pi/piNpiS/indi_fasta/${indi}.fasta -o /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/effective_population_size/pi/piNpiS/indi_fasta/${indi}.rename.fasta
python 1_extract_CDS_from_chromsome.py --individual ${indi}
python 2_Obtain_diploid_sequences.py --individual ${indi}