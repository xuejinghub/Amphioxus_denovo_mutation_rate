#!/bin/bash
#SBATCH --job-name=platanus
#SBATCH -c 20
#SBATCH --time=24:00:00                     
#SBATCH --mem=100G

prefix=Bf-${SLURM_ARRAY_TASK_ID}                          

mkdir whole_genome_assembly/${prefix}

cd whole_genome_assembly/${prefix}

mkdir tmp

fastq_dictory=fastq

# # Step 1: Filter reads with Fastp, removing reads with >10% N and minimum quality of 20
if [ ! -e $fastq_dictory/${prefix}_qc1_1.fq.gz ] || [ ! -e $fastq_dictory/${prefix}_qc1_2.fq.gz]
then
    fastp -i $fastq_dictory/${prefix}_1.fq.gz -I $fastq_dictory/${prefix}_2.fq.gz -o $fastq_dictory/${prefix}_qc1_1.fq.gz \
    -O $fastq_dictory/${prefix}_qc1_2.fq.gz \
      --n_base_limit 5 \
      --length_required 15 \
      -q 20 \
      -u 40 \
      -c \
      --detect_adapter_for_pe \
      --adapter_sequence=AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
      --adapter_sequence_r2=AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
      --failed_out tmp/${prefix}_qc1_failed_reads.fastq \
      -j tmp/${prefix}_qc1_fastp_report.json \
      -h tmp/${prefix}_qc1_fastp_report.html \
      -w 4
fi 

# # Step 2: Filter low-quality reads with Fastp
if [ ! -e $fastq_dictory/${prefix}_qc2_1.fq.gz ] || [ ! -e $fastq_dictory/${prefix}_qc2_2.fq.gz]
then
    fastp -i $fastq_dictory/${prefix}_qc1_1.fq.gz -I $fastq_dictory/${prefix}_qc1_2.fq.gz \
        -o $fastq_dictory/${prefix}_qc2_1.fq.gz -O $fastq_dictory/${prefix}_qc2_2.fq.gz \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --failed_out tmp/${prefix}_qc2_failed_reads.fastq \
        --detect_adapter_for_pe \
        -j tmp/${prefix}_qc2_fastp_report.json \
        -h tmp/${prefix}_qc2_fastp_report.html \
        -c
fi
# # Step 3: Assemble reads with Platanus

platanus_allee assemble \
    -m 200 -t 20 \
    -tmp tmp \
    -o ${prefix}_platanus \
    -f $fastq_dictory/${prefix}_qc2_1.fq.gz $fastq_dictory/${prefix}_qc2_2.fq.gz -o ${prefix}.platanus

platanus_allee phase \
    -t 20 \
    -i 3 \
    -c ${prefix}.platanus_contig.fa ${prefix}.platanus_32merFrq.tsv \
    -o ${prefix}_platanus \
    -IP1 $fastq_dictory/${prefix}_qc2_1.fq.gz $fastq_dictory/${prefix}_qc2_2.fq.gz


seqkit seq -m 150 ${prefix}_platanus_allPhasedScaffold.fa > ${prefix}_platanus_allPhasedScaffold.min150.fa

quast -m 150 -o ${prefix}_platanus_quast -t 20 ${prefix}_platanus_allPhasedScaffold.min150.fa
