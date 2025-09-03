#!/bin/bash
#SBATCH -w node80
#SBATCH -c 3
#SBATCH --time=48:00:00                     
#SBATCH --mem=10G
#SBATCH -J CombineGVCFs
#SBATCH -o log/CombineGVCFs_%a.%j.out

if [ ! -d gvcf ]; then
    mkdir gvcf
fi

prefix=Chr${SLURM_ARRAY_TASK_ID}

if [ "$prefix" = "Chr16" ]; then
    prefix="Chr1620"
fi

gatk --java-options "-Xmx10g" CombineGVCFs \
    -R bf.dnm/reference/bf.fa \
    --tmp-dir bf.snvCalling/combinedgVCF/tmp \
    -O gvcf/Bf_inds35_${prefix}.g.vcf.gz \
    -L ${prefix} \
	-V haplotypecaller/SRR12010215/SRR12010215.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010217/SRR12010217.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010218/SRR12010218.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010219/SRR12010219.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010220/SRR12010220.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010222/SRR12010222.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010224/SRR12010224.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010226/SRR12010226.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010229/SRR12010229.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010230/SRR12010230.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010231/SRR12010231.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010233/SRR12010233.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010234/SRR12010234.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010235/SRR12010235.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010237/SRR12010237.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010239/SRR12010239.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010240/SRR12010240.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010241/SRR12010241.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010242/SRR12010242.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010243/SRR12010243.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010244/SRR12010244.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010247/SRR12010247.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010252/SRR12010252.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010253/SRR12010253.haplotypecaller.g.vcf.gz \
	-V haplotypecaller/SRR12010254/SRR12010254.haplotypecaller.g.vcf.gz \
	-V variant_calling/haplotypecaller/SRR12010221/SRR12010221.haplotypecaller.g.vcf.gz \
	-V variant_calling/haplotypecaller/SRR12010225/SRR12010225.haplotypecaller.g.vcf.gz \
	-V variant_calling/haplotypecaller/SRR12010228/SRR12010228.haplotypecaller.g.vcf.gz \
	-V variant_calling/haplotypecaller/SRR12010232/SRR12010232.haplotypecaller.g.vcf.gz \
	-V variant_calling/haplotypecaller/SRR12010236/SRR12010236.haplotypecaller.g.vcf.gz \
	-V variant_calling/haplotypecaller/SRR12010246/SRR12010246.haplotypecaller.g.vcf.gz \
	-V variant_calling/haplotypecaller/SRR12010248/SRR12010248.haplotypecaller.g.vcf.gz \
	-V variant_calling/haplotypecaller/SRR12010251/SRR12010251.haplotypecaller.g.vcf.gz \
	-V variant_calling/haplotypecaller/SRR12010250/SRR12010250.haplotypecaller.g.vcf.gz \
	-V variant_calling/haplotypecaller/SRR12010223/SRR12010223.haplotypecaller.g.vcf.gz


gatk --java-options "-Xmx20000M -XX:ParallelGCThreads=3" GenotypeGVCFs \
  --variant gvcf/Bf_inds35_${prefix}.g.vcf.gz \
  --heterozygosity 0.3 \
  --indel-heterozygosity 0.03 \
  -L ${prefix} \
  --include-non-variant-sites \
  --output chrm_vcf/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.vcf.gz \
  -R bf.dnm/reference/bf.fa \
  --tmp-dir tmp 

gatk --java-options "-Xmx20000M -XX:ParallelGCThreads=3" SelectVariants \
    -V chrm_vcf/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.vcf.gz \
    -select-type SNP \
    -O chrm_vcf/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.SNPs.vcf.gz

gatk --java-options "-Xmx20000M -XX:ParallelGCThreads=3" VariantFiltration \
    --variant chrm_vcf/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.SNPs.vcf.gz  \
    --filter "QD < 2.0" --filter-name "QD2" \
    --filter "QUAL < 30.0" --filter-name "QUAL30" \
    --filter "SOR > 3.0" --filter-name "SOR3" \
    --filter "FS > 60.0" --filter-name "FS60" \
    --filter "MQ < 40.0" --filter-name "MQ40" \
    --filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --output chrm_vcf/hardfilter/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.SNPs.filtering.vcf.gz

bcftools view \
    -f PASS \
    --threads 3 \
    -O z -o chrm_vcf/hardfilter/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.SNPs.hardfiltered.vcf.gz \
    chrm_vcf/hardfilter/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.SNPs.filtering.vcf.gz

gatk --java-options "-Xmx20000M -XX:ParallelGCThreads=3" SelectVariants \
    -V chrm_vcf/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.vcf.gz \
    -select-type INDEL \
    -O chrm_vcf/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.INDELs.vcf.gz

gatk --java-options "-Xmx20000M -XX:ParallelGCThreads=3" VariantFiltration \
    --variant chrm_vcf/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.INDELs.vcf.gz  \
    --filter "QD < 2.0" --filter-name "QD2" \
    --filter "QUAL < 30.0" --filter-name "QUAL30" \
    --filter "SOR > 10.0" --filter-name "SOR10" \
    --filter "FS > 200.0" --filter-name "FS200" \
    --filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    --output chrm_vcf/hardfilter/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.INDELs.filtering.vcf.gz 

bcftools view \
    -f PASS \
    --threads 3 \
    -O z -o chrm_vcf/hardfilter/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.INDELs.hardfiltered.vcf.gz \
    chrm_vcf/hardfilter/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.INDELs.filtering.vcf.gz

# 新增非变异位点过滤
gatk --java-options "-Xmx20000M -XX:ParallelGCThreads=3" SelectVariants \
    -V chrm_vcf/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.vcf.gz \
    --select-type-to-include NO_VARIATION \
    -O chrm_vcf/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.INVARIANT.vcf.gz

bcftools +fill-tags chrm_vcf/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.INVARIANT.vcf.gz  -t F_MISSING | bcftools view -i 'F_MISSING<0.25' -Oz -o chrm_vcf/hardfilter/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.INVARIANT.hardfiltered.vcf.gz

tabix -p vcf chrm_vcf/hardfilter/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.INVARIANT.hardfiltered.vcf.gz

tabix -p vcf chrm_vcf/hardfilter/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.SNPs.hardfiltered.vcf.gz

bcftools concat \
    --allow-overlaps \
    chrm_vcf/hardfilter/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.SNPs.hardfiltered.vcf.gz \
    chrm_vcf/hardfilter/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_variants.INVARIANT.hardfiltered.vcf.gz \
    | bcftools sort - \
    | bcftools view -Oz \
    -o chrm_vcf/hardfilter/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_allsites.hardfiltered.vcf.gz \
    -

tabix -p vcf chrm_vcf/hardfilter/Amphioxus_bf_35inds_autosome_highdepth_${prefix}_allsites.hardfiltered.vcf.gz


