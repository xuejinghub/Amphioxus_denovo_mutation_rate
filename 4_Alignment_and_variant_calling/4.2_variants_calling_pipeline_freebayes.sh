#!/bin/bash
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=freebayes
#SBATCH --output=log/freebayes_%j_%a.out

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

id=${SLURM_ARRAY_TASK_ID}
prefix=Bf-${id}
INDEX=reference/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa
startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`

source /public/home/xuejing/mambaforge/etc/profile.d/conda.sh
conda activate freebayes

# Step 0: Check if the bed of CG is generated
if [ ! -e bed/Bf-${id}_platanus_i3_allPhasedScaffold.rename.min150bp.fa.bed ]
then
    awk -F "\t" 'NR==FNR {a[$1]=1;next} $1 in a{print $0}' callable_genome/tmp/Bf-${id}.contig_inheritance_v0.8.postprocess.list callable_genome/progeny_depth/result/Bf-${id}.contig_inheritance_v0.8.postprocess.remove3foldsites.sort.bed > bed/Bf-${id}_platanus_i3_allPhasedScaffold.rename.min150bp.fa.bed
    mkdir tmp/Bf-${id}
    split -d --numeric-suffixes=10 -l 2200 tmp/Bf-${id}_platanus_i3_allPhasedScaffold.rename.min150bp.fa.bed tmp/Bf-${id}/Bf-${id}_platanus_i3_allPhasedScaffold.rename.min150bp.fa.bed.part
fi

## Step 1: Freebayes calling variants
freebayes -f ${INDEX} \
    -m 0 \
    --legacy-gls \
    -t callable_genome/bed/${prefix}.contig_inheritance_v0.8.postprocess.remove_similarity.remove3foldsites.sorted.bed \
    bam/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam \
    | bgzip -c > vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz
tabix -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz

source /public/home/xuejing/mambaforge/etc/profile.d/conda.sh
conda activate base

## Step 2: Filter variants
bcftools isec \
    -n~100 \
    -c all \
    -p tmp/${prefix}.out \
    vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz  \
    vcf/Bf_M.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz \
    vcf/Bf_P.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz


cut -f1,2 tmp/${prefix}.out/sites.txt > tmp/${prefix}.out/sites.list
bcftools view  \
    -R tmp/${prefix}.out/sites.list \
    vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz -Ov | \
    bcftools norm -m -any | \
    bcftools norm -a -f ${INDEX} | \
    bcftools view  \
    -v snps -e '(INFO/DP)-(INFO/AO)>2 || INFO/DP <= 10' \
    -Oz -o vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.vcf.gz - \
    && tabix -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.vcf.gz
python script/filter_vcf.py --in vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.vcf.gz --bam ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/Bf_P.mem2..md.uniq.bam --bam ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/Bf_M.mem2..md.uniq.bam --out vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.filtered.vcf

cat vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.filtered.vcf | grep -v "#" |  grep "TYPE=snp" > result/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.filtered.vcf


source /public/home/xuejing/mambaforge/etc/profile.d/conda.sh
conda activate freebayes

## Step 3: Freebayes calling unphased variants 
if [ ! -e vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.unphased.vcf.gz.tbi ]
then
freebayes -f ${INDEX} \
    -m 0 \
    --legacy-gls \
    -t ${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/similarity_analysis/postprocess_result/${prefix}/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.uncallableE3.bed \
    ${user_home}/Amphioxus/new_align/reassembly/read_align/default/bam/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam \
    | bgzip -c > vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.unphased.vcf.gz
tabix -f -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.unphased.vcf.gz
fi

source /public/home/xuejing/mambaforge/etc/profile.d/conda.sh
conda activate base

## Step 4: Filter unphased variants
bcftools isec \
    -n~100 \
    -c all \
    -p tmp/${prefix}.unphased.out \
    vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.unphased.vcf.gz  \
    vcf/Bf_M.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz \
    vcf/Bf_P.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz

cut -f1,2 tmp/${prefix}.unphased.out/sites.txt > tmp/${prefix}.unphased.out/sites.list
bcftools view  \
    -R tmp/${prefix}.unphased.out/sites.list \
    vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.unphased.vcf.gz -Ov | \
    bcftools norm -m -any | \
    bcftools norm -a -f ${INDEX} | \
    bcftools view  \
    -v snps -e '(INFO/DP)-(INFO/AO)>2 || INFO/DP <= 10' \
    -Oz -o vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.vcf.gz - \
    && tabix -f -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.vcf.gz
python script/filter_vcf.py --in vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.vcf.gz --bam ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/Bf_P.mem2..md.uniq.bam --bam ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/Bf_M.mem2..md.uniq.bam --out vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.filtered.vcf

cat vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.filtered.vcf | grep -v "#" |  awk  '$0!~/TYPE=indel/' > result/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.filtered.vcf

endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
