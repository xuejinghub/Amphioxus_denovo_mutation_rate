#!/bin/bash
#SBATCH -w node80
#SBATCH -c 20
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=QC_preprocess_mpileup
#SBATCH --output=log/QC_preprocess_%a.out

if [ ! -d tmp ]
then
    mkdir tmp
fi
if [ ! -d bam ]
then    
    mkdir bam
fi
if [ ! -d vcf ]
then
    mkdir vcf
fi
if [ ! -d report ]
then
    mkdir report
fi
if [ ! -d result ]
then
    mkdir result
fi
count_mutation=${1}
prefix=Bf${SLURM_ARRAY_TASK_ID}_mut${1}
fastq_dictory=FNR/simulation_data
align_software=mem2
parameters="Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp"
threads=20
INDEX=reference/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa
if [ ! -d tmp/${prefix}.out ]
then
    mkdir tmp/${prefix}.out
fi
if [ ! -e ${INDEX}.amb ]
then    
    bwa-mem2 index $INDEX
fi

# # Step 1: Filter reads with Fastp, removing reads with >10% N and minimum quality of 20
if [ ! -e $fastq_dictory/${prefix}_qc1_1.fq.gz ] || [ ! -e $fastq_dictory/${prefix}_qc1_2.fq.gz ]
then
    fastp -i $fastq_dictory/${prefix}1.fq.gz -I $fastq_dictory/${prefix}2.fq.gz -o $fastq_dictory/${prefix}_qc1_1.fq.gz \
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
if [ ! -e $fastq_dictory/${prefix}_qc2_1.fq.gz ] || [ ! -e $fastq_dictory/${prefix}_qc2_2.fq.gz ] 
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

if [ ! -d report/${prefix} ]
then
    mkdir -p report/${prefix}
fi

# Step 3: Align filtered reads with BWA
if [ ! -e bam/${prefix}.${align_software}.${parameters}.bam ]
then
    if [ "$align_software" == "mem2" ]; then
        bwa-mem2 mem \
            -R "@RG\tID:null.${prefix}.1\tPU:1\tSM:${prefix}\tLB:${prefix}\tDS:$INDEX\tPL:ILLUMINA" \
            -t ${threads} \
            -M -k 19 \
            $INDEX \
            $fastq_dictory/${prefix}_qc2_1.fq.gz $fastq_dictory/${prefix}_qc2_2.fq.gz |  \
            samtools sort -@ ${threads} -o bam/${prefix}.${align_software}.${parameters}.bam -
    elif [ "$align_software" == "bwa" ]; then
        bwa mem \
            -R "@RG\tID:null.${prefix}.1\tPU:1\tSM:${prefix}\tLB:${prefix}\tDS:$INDEX\tPL:ILLUMINA" \
            -t ${threads} \
            -M -k 19 \
            $INDEX \
            $fastq_dictory/${prefix}_qc2_1.fq.gz $fastq_dictory/${prefix}_qc2_2.fq.gz | \
            samtools sort -@ ${threads} -o bam/${prefix}.${align_software}.${parameters}.bam -
    fi
fi

if [ ! -e bam/${prefix}.${align_software}.${parameters}.bam.bai ]
then
    samtools index bam/${prefix}.${align_software}.${parameters}.bam
fi

# Step 4: Mark duplicates using Picard
if [ ! -e bam/${prefix}.${align_software}.${parameters}.md.bam ]
then
    gatk MarkDuplicates \
        I=bam/${prefix}.${align_software}.${parameters}.bam \
        O=bam/${prefix}.${align_software}.${parameters}.md.bam \
        M=bam/${prefix}.${align_software}.${parameters}.md.metrics.txt \
        REMOVE_DUPLICATES=true # Set to true if you want to remove duplicates, or false to keep them marked

    samtools index bam/${prefix}.${align_software}.${parameters}.md.bam
fi

# Step 5: Filter unmapped reads
if [ ! -e bam/${prefix}.${align_software}.${parameters}.md.uniq.bam ]
then
    sambamba view -t ${threads} -h -f bam \
        -F "not (secondary_alignment or supplementary)"   \
        -p -l 9 \
        bam/${prefix}.${align_software}.${parameters}.md.bam \
        -o bam/${prefix}.${align_software}.${parameters}.md.uniq.bam
fi

if [ ! -e bam/${prefix}.${align_software}.${parameters}.md.uniq.bam.bai ]
then
    samtools index bam/${prefix}.${align_software}.${parameters}.md.uniq.bam
fi

# Step 6: Index
if [ ! -e bam/${prefix}.${align_software}.${parameters}.md.uniq.bam.depth.mosdepth.summary.txt ]
then
    mosdepth bam/${prefix}.${align_software}.${parameters}.md.uniq.bam.depth bam/${prefix}.${align_software}.${parameters}.md.uniq.bam
fi



