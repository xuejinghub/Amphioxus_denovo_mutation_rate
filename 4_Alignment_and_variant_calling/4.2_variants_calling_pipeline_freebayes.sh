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

## Step 2: Filter variants
bcftools isec \
    -n~100 \
    -c all \
    -p tmp/${prefix}.out vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz  \
    Bf_P.mem2.default_min150_phased.md.uniq.bcftools.vcf.gz \
    Bf_M.mem2.default_min150_phased.md.uniq.bcftools.vcf.gz

cut -f1,2 tmp/${prefix}.out/sites.txt > tmp/${prefix}.out/sites.list
bcftools view  \
    -R tmp/${prefix}.out/sites.list \
    vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz -Ov |\
    bcftools view  \
    -v snps -e '(INFO/DP)-(INFO/AO)>2 || INFO/DP <= 10' \
    -Oz -o vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.vcf.gz - \
    && tabix -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.vcf.gz
bcftools norm -f ${INDEX} vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.vcf.gz -Oz -o vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.vcf.gz && tabix -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.vcf.gz
python filter_vcf.py --in vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.vcf.gz --bam bam/Bf_P.mem2..md.uniq.bam --bam bam/Bf_M.mem2..md.uniq.bam --out vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.filtered.vcf

cat vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.filtered.vcf | grep -v "#" |  grep "TYPE=snp" > result/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.filtered.vcf


endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
