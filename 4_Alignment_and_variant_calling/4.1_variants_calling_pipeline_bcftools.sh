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

prefix=Bf-${SLURM_ARRAY_TASK_ID}
fastq_dictory=fastq
align_software=mem2
parameters="Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp"
threads=20
INDEX=Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa
if [ ! -e ${INDEX}.amb ]
then    
    bwa-mem2 index $INDEX
fi
## Step 0: Check if the input files exist
if [ ! -e $fastq_dictory/${prefix}_1.fq.gz ] || [ ! -e $fastq_dictory/${prefix}_2.fq.gz ]
then
    fastq1=`readlink -f bf.data/Bf-${SLURM_ARRAY_TASK_ID}/*/*_1.fq.gz`
    fastq2=`readlink -f bf.data/Bf-${SLURM_ARRAY_TASK_ID}/*/*_2.fq.gz`
    ln -s $fastq1 $fastq_dictory/${prefix}_1.fq.gz
    ln -s $fastq2 $fastq_dictory/${prefix}_2.fq.gz
fi

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

if [ ! -d report/${prefix} ]
then
    mkdir -p report/${prefix}
fi

if [ ! -e report/${prefix}/${prefix}_qc1_fastp_report.html ]
then
    fastqc -o report/${prefix} $fastq_dictory/${prefix}_qc1_1.fq.gz $fastq_dictory/${prefix}_qc1_2.fq.gz -t ${threads}
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

# Step 7: Call variants with BCFtools
if [ ! -e vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.bcftools.vcf.gz ]
then
    bcftools mpileup \
        -a FORMAT/SP,FORMAT/ADF,FORMAT/ADR,INFO/SCR,FORMAT/QS,FORMAT/AD,FORMAT/DP,FORMAT/SCR,FORMAT/NMBZ \
        -f ${INDEX} \
        --threads ${threads} \
        -Ou bam/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam | bcftools call -mv -Oz \
        -o vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.bcftools.vcf.gz
    tabix -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.bcftools.vcf.gz
fi

if [ ! -e vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.vcf.gz.tbi ]
then
    bcftools isec -n~100 -c all -p ${prefix}.out vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.bcftools.vcf.gz vcf/Bf_P.mem2.default_min150_phased.md.uniq.bcftools.vcf.gz vcf/Bf_M.mem2.default_min150_phased.md.uniq.bcftools.vcf.gz
    cut -f1,2 ${prefix}.out/sites.txt > ${prefix}.out/sites.list
    bcftools view -R ${prefix}.out/sites.list -e 'GT="het" || INFO/DP <= 10' -Oz -o vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.vcf.gz vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.bcftools.vcf.gz
    tabix -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.vcf.gz
fi
if [ ! -e vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.filtered.vcf ]
then
    python filter_vcf.py --in vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.vcf.gz --bam reassembly/read_align/reference/read_align/bam/Bf_P.mem2..md.uniq.bam --bam reassembly/read_align/reference/read_align/bam/Bf_M.mem2..md.uniq.bam --out vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.filtered.vcf
fi

cat vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.filtered.vcf | grep -v "#" |  grep -v "INDEL" | awk -F "\t" '$6>225' > result/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.filtered.vcf
mawk -F "\t" 'NR==FNR { 
    a[$0] = 1
    next
} 
    $1 in a { 
    print 
}' inheritance_filter/pass_list/${prefix}.passlist result/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.filtered.vcf > result2/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.filtered.pass.vcf

