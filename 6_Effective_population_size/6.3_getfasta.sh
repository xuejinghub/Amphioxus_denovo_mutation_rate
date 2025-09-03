#!/bin/bash
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=getfasta
#SBATCH --output=log/getfasta_%j.out

prefix=Chr${SLURM_ARRAY_TASK_ID}

if [ "$prefix" = "Chr16" ]; then
    prefix="Chr1620"
fi

if [ ! -d vcf ]
then
    mkdir vcf
fi
indi=$(head -n ${SLURM_ARRAY_TASK_ID} indi.list | tail -1)
gatk FastaAlternateReferenceMaker -R reference/bf.fa -O indi_fasta/${indi}.fasta --use-iupac-sample ${indi} -V effective_population_size/pi/piNpiS/Amphioxus_bf_35inds_autosome_highdepth_allchrom_variantSites.biallele.hardfiltered.vcf.gz
