#!/bin/bash
#SBATCH -w node80
#SBATCH -c 20
#SBATCH --time=48:00:00                     
#SBATCH --mem=100G
#SBATCH --job-name=sim80
#SBATCH --output=log/sim80_%j.out

chr=${1}
het=$3
parent=$2
depth=40
species=bf
prefix=${chr}_sim_het${het}_${parent}_dep${depth}_batch${4}
genome_size=$(cat ${user_home}/Amphioxus/bf.dnm/reference/bf.fa.fai | awk -v chr="${chr}" '$1==chr{print $2}')
refrence=${user_home}/Amphioxus/bf.dnm/reference/bf.fa
cpu=20
fastq_dictory=wd/${prefix}/fq
min_score=200
align_software=mem2
parameters=${chr}_sim_het${het}_${parent}_dep${depth}_allPhasedScaffold.min150bp
threads=${cpu}

simg_wd=${chr}_sim_het${het}_${parent}_dep${depth}_batch${4}
sim_genome_dir=sim_genome/${simg_wd}
mkdir ${sim_genome_dir}
cd sim_genome

bcftools view -r ${chr} -Oz bf.snp.vcf.gz -o ${simg_wd}/bf.snp.${chr}.vcf.gz

tabix -p vcf ${simg_wd}/bf.snp.${chr}.vcf.gz

python simulate_inhom_poisson_v0.4.py \
    --ref /public5/home/taolei/bf/review/supp/r1.1/${chr}.fa \
    --win 1000 \
    --obs-vcf ${simg_wd}/bf.snp.${chr}.vcf.gz \
    --heterozygosity 0.03 \
    --poisson \
    --hap1 ${simg_wd}/${chr}_het${het}_simv0.4_hap1.fa --hap2 ${simg_wd}/${chr}_het${het}_simv0.4_hap2.fa \
    --vcf ${simg_wd}/genome_simulated_variants.vcf.gz \
    --seed 199${4} \
    --plot-prefix ${simg_wd}/sim \
    --nonvar-total 3000 \
    --nonvar-bed ${simg_wd}/nonvariant_regions.bed 
    # --nonvar-mutate-frac 0.1 \
    # --mut-nonvar-bed ${wd}/mutated_nonvariant_regions.bed \
    # --mut-nonvar-variants-bed ${wd}/mutated_nonvariant_variants.bed \
    # --out-count-dist $wd/win_count_distribution.tsv

cd ..

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

sed "s/${chr}/${chr}_hap1/g" ${sim_genome_dir}/${chr}_het${het}_simv0.4_hap1.fa > ${sim_genome_dir}/${chr}_sim_het${het}_${parent}_dep${depth}_hap1.fa
sed "s/${chr}/${chr}_hap2/g" ${sim_genome_dir}/${chr}_het${het}_simv0.4_hap2.fa > ${sim_genome_dir}/${chr}_sim_het${het}_${parent}_dep${depth}_hap2.fa

cat ${sim_genome_dir}/${chr}_sim_het${het}_${parent}_dep${depth}_hap1.fa ${sim_genome_dir}/${chr}_sim_het${het}_${parent}_dep${depth}_hap2.fa > fa/${chr}_sim_het${het}_${parent}_dep${depth}_poisson_diploid.fa

if [ ! -e $fastq_dictory/${prefix}_1.fq.gz ] || [ ! -e $fastq_dictory/${prefix}_2.fq.gz ]
then
art_illumina -p -m 250 -s 50 -f $depth -l 150 -ef -ss HS25 -i fa/${chr}_sim_het${het}_${parent}_dep${depth}_poisson_diploid.fa -o fq/${prefix}_
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

cd fa

${user_home}/Amphioxus/new_align/reassembly/platanus_allee assemble \
  -m 200 -t ${cpu} -tmp ../tmp \
  -f $fastq_dictory/${prefix}_qc2_1.fq.gz $fastq_dictory/${prefix}_qc2_2.fq.gz \
  -o ${chr}_sim_het${het}_${parent}_dep${depth}
${user_home}/Amphioxus/new_align/reassembly/platanus_allee phase \
  -t ${cpu} -i 3 \
  -c ${chr}_sim_het${het}_${parent}_dep${depth}_contig.fa ${chr}_sim_het${het}_${parent}_dep${depth}_32merFrq.tsv \
  -IP1 $fastq_dictory/${prefix}_qc2_1.fq.gz $fastq_dictory/${prefix}_qc2_2.fq.gz \
  -o ${chr}_sim_het${het}_${parent}_dep${depth}

quast ${chr}_sim_het${het}_${parent}_dep${depth}_allPhasedScaffold.fa

cd ..

if [ ! -e fa/${chr}_sim_het${het}_${parent}_dep${depth}_allPhasedScaffold.min150.fa ]
then
seqkit seq -m 150 fa/${chr}_sim_het${het}_${parent}_dep${depth}_allPhasedScaffold.fa > fa/${chr}_sim_het${het}_${parent}_dep${depth}_allPhasedScaffold.min150.fa
fi

INDEX=fa/${chr}_sim_het${het}_${parent}_dep${depth}_allPhasedScaffold.min150.fa

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
