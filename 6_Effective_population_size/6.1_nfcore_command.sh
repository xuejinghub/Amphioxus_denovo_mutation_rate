#!/bin/bash
#SBATCH -w node80
#SBATCH -c 50
#SBATCH --time=48:00:00                     
#SBATCH --mem=400G
#SBATCH --job-name=sarek
#SBATCH --output=log/sarek_%j.out

nextflow run sarek -profile singularity --input bf_input.batch1_1.csv --fasta reference/bf.fa --igenomes_ignore --genome custom --max_memory 400.GB --tools haplotypecaller  --outdir  batch1_bf.out --max_cpus 50  --aligner bwa-mem2 --skip_tools baserecalibrator -resume disturbed_kalam --joint_germline -c sarek/conf/modules/custom.config
