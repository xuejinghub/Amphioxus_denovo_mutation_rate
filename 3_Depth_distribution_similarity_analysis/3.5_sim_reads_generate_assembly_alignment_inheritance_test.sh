#!/bin/bash
#SBATCH -w node80
#SBATCH -c 20
#SBATCH --time=48:00:00                     
#SBATCH --mem=100G
#SBATCH --job-name=1sim_inheritance
#SBATCH --output=log/1sim_inheritance_%j.out

chr=Chr3
het=0.03
depth=40
species=bf
batch=${1}
prefix=${chr}_sim_het${het}_F1batch${batch}_dep${depth}
genome_size=$(cat ${user_home}/Amphioxus/bf.dnm/reference/bf.fa.fai | awk -v chr="${chr}" '$1==chr{print $2}')
refrence=${user_home}/Amphioxus/bf.dnm/reference/bf.fa
cpu=20
fastq_dictory=${user_home}/Amphioxus/new_align/reassembly/read_align/cosin_similarity/simulation/wd/${prefix}/fq
min_score=200
align_software=mem2
parameters=${chr}_sim_het${het}_MP_dep${depth}_allPhasedScaffold.min150.rename
threads=${cpu}

sim_genome_dir=${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/align_readdepth/sim/pisson_assebly/sim_genome

if [ ! -d wd/${prefix} ]; then
    mkdir -p wd/${prefix}
fi

cd wd/${prefix}

if [ ! -d tmp ]; then
    mkdir tmp
fi

if [ ! -d fq ]; then
    mkdir fq
fi

if [ ! -d fa ]; then
    mkdir fa
fi

if [ ! -d bam ]; then
    mkdir bam
fi

if [ ! -d align ]; then
    mkdir align
fi

if [ ! -d result ]; then
    mkdir result
fi

source ${conda_home}/etc/profile.d/conda.sh
conda activate art

sed "s/${chr}/${chr}_P_hap${batch}/g" ${sim_genome_dir}/${chr}_sim_het${het}_P_dep${depth}_batch1/${chr}_het${het}_simv0.4_hap${batch}.fa > fa/${chr}_sim_het${het}_P_dep${depth}_hap${batch}.fa
sed "s/${chr}/${chr}_M_hap${batch}/g" ${sim_genome_dir}/${chr}_sim_het${het}_M_dep${depth}_batch2/${chr}_het${het}_simv0.4_hap${batch}.fa > fa/${chr}_sim_het${het}_M_dep${depth}_hap${batch}.fa
cat fa/${chr}_sim_het${het}_P_dep${depth}_hap${batch}.fa fa/${chr}_sim_het${het}_M_dep${depth}_hap${batch}.fa > fa/${chr}_sim_het${het}_F1_dep${depth}_poisson_diploid.fa

if [ ! -e $fastq_dictory/${prefix}_1.fq.gz ] || [ ! -e $fastq_dictory/${prefix}_2.fq.gz ]
then
art_illumina -p -m 250 -s 50 -f $depth -l 150 -ef -ss HS25 -i fa/${chr}_sim_het${het}_F1_dep${depth}_poisson_diploid.fa -o fq/${prefix}_
fi

source ${conda_home}/etc/profile.d/conda.sh
conda activate align

if [ ! -e $fastq_dictory/${prefix}_1.fq.gz ] || [ ! -e $fastq_dictory/${prefix}_2.fq.gz ]
then
    bgzip -c $fastq_dictory/${prefix}_1.fq > $fastq_dictory/${prefix}_1.fq.gz
    bgzip -c $fastq_dictory/${prefix}_2.fq > $fastq_dictory/${prefix}_2.fq.gz
fi

# # Step 1: Filter reads with Fastp, removing reads with >10% N and minimum quality of 20
if [ ! -e $fastq_dictory/${prefix}_qc1_1.fq.gz ] || [ ! -e $fastq_dictory/${prefix}_qc1_2.fq.gz ]
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


if [ ! -e fa/${chr}_sim_het${het}_MP_dep${depth}_allPhasedScaffold.min150.rename.fa ]
then
cat  ${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/align_readdepth/sim/pisson_assebly/wd/${chr}_sim_het${het}_P_dep40_batch1/fa/${chr}_sim_het${het}_P_dep40_allPhasedScaffold.min150.fa | awk '$0~/^>/{print $0"_Bf_P"} $0!~/^>/ {print;}' > fa/${chr}_sim_het${het}_P_dep40_allPhasedScaffold.min150.rename.fa
cat  ${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/align_readdepth/sim/pisson_assebly/wd/${chr}_sim_het${het}_M_dep40_batch2/fa/${chr}_sim_het${het}_M_dep40_allPhasedScaffold.min150.fa | awk '$0~/^>/{print $0"_Bf_M"} $0!~/^>/ {print;}' > fa/${chr}_sim_het${het}_M_dep40_allPhasedScaffold.min150.rename.fa
cat fa/${chr}_sim_het${het}_P_dep40_allPhasedScaffold.min150.rename.fa fa/${chr}_sim_het${het}_M_dep40_allPhasedScaffold.min150.rename.fa > fa/${chr}_sim_het${het}_MP_dep${depth}_allPhasedScaffold.min150.rename.fa
fi

INDEX=fa/${chr}_sim_het${het}_MP_dep${depth}_allPhasedScaffold.min150.rename.fa

if [ ! -e ${INDEX}.amb ]
then    
    bwa-mem2 index $INDEX
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
            $fastq_dictory/${prefix}_qc2_1.fq.gz $fastq_dictory/${prefix}_qc2_2.fq.gz | \
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

#######################maternal align########################################

fastq_dictory=${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/align_readdepth/sim/pisson_assebly/wd/Chr3_sim_het0.03_M_dep40_batch2/fq
prefix="Chr3_sim_het0.03_M_dep40_batch2"

# Step 3: Align filtered reads with BWA
if [ ! -e bam/${prefix}.${align_software}.${parameters}.bam ]
then
    if [ "$align_software" == "mem2" ]; then
        bwa-mem2 mem \
            -R "@RG\tID:null.${prefix}.1\tPU:1\tSM:${prefix}\tLB:${prefix}\tDS:$INDEX\tPL:ILLUMINA" \
            -t ${threads} \
            -M -k 19 \
            $INDEX \
            $fastq_dictory/${prefix}_qc2_1.fq.gz $fastq_dictory/${prefix}_qc2_2.fq.gz | \
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

#######################paternal align########################################

fastq_dictory=${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/align_readdepth/sim/pisson_assebly/wd/Chr3_sim_het0.03_P_dep40_batch1/fq/
prefix="Chr3_sim_het0.03_P_dep40_batch1"

# Step 3: Align filtered reads with BWA
if [ ! -e bam/${prefix}.${align_software}.${parameters}.bam ]
then
    if [ "$align_software" == "mem2" ]; then
        bwa-mem2 mem \
            -R "@RG\tID:null.${prefix}.1\tPU:1\tSM:${prefix}\tLB:${prefix}\tDS:$INDEX\tPL:ILLUMINA" \
            -t ${threads} \
            -M -k 19 \
            $INDEX \
            $fastq_dictory/${prefix}_qc2_1.fq.gz $fastq_dictory/${prefix}_qc2_2.fq.gz | \
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


