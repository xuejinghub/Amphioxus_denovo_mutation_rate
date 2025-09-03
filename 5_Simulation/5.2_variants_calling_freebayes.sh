#!/bin/bash
#SBATCH -w node80
#SBATCH -c 2
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=freebayes
#SBATCH --output=log/freebayes_%j_%a.out
id=$1
# count_mutation=${1}
prefix=Bf${id}_mut${2}
# fastq_dictory=FNR/simulation_data
fastq_dictory=FNR/simulation_data
align_software=mem2
parameters="Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp"
threads=50
INDEX=reference/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa
startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`


if [ ! -e tmp/${prefix}/${prefix}_platanus_i3_allPhasedScaffold.rename.min150bp.fa.bed.part10 ]
then
    awk -F "\t" 'NR==FNR {a[$1]=1;next} $1 in a{print $0}' FNR/callable_genome/tmp/Bf-${id}.contig_inheritance_v0.8.postprocess.list  default/callable_genome/progeny_depth/result/Bf-${id}.contig_inheritance_v0.8.postprocess.remove3foldsites.sort.bed > tmp/${prefix}_platanus_i3_allPhasedScaffold.rename.min150bp.fa.bed
    mkdir tmp/${prefix}
    split -d --numeric-suffixes=10 -l 2200 tmp/${prefix}_platanus_i3_allPhasedScaffold.rename.min150bp.fa.bed tmp/${prefix}/${prefix}_platanus_i3_allPhasedScaffold.rename.min150bp.fa.bed.part
fi

freebayes -f ${INDEX} \
    -m 0 \
    --legacy-gls \
    -t tmp/${prefix}/${prefix}_platanus_i3_allPhasedScaffold.rename.min150bp.fa.bed.part${SLURM_ARRAY_TASK_ID} \
    FNR/bam/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam \
    | bgzip -c > vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.part${SLURM_ARRAY_TASK_ID}.vcf.gz && tabix -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.part${SLURM_ARRAY_TASK_ID}.vcf.gz
bcftools isec \
    -n~100 \
    -c all \
    -p tmp/${prefix}.out/part${SLURM_ARRAY_TASK_ID} vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.part${SLURM_ARRAY_TASK_ID}.vcf.gz \
    default/vcf/Bf_P.mem2.default_min150_phased.md.uniq.bcftools.vcf.gz \
    default/vcf/Bf_M.mem2.default_min150_phased.md.uniq.bcftools.vcf.gz
cut -f1,2 tmp/${prefix}.out/part${SLURM_ARRAY_TASK_ID}/sites.txt > tmp/${prefix}.out/part${SLURM_ARRAY_TASK_ID}/sites.list
bcftools view -R tmp/${prefix}.out/part${SLURM_ARRAY_TASK_ID}/sites.list vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.part${SLURM_ARRAY_TASK_ID}.vcf.gz | bcftools view  -e '(INFO/DP)-(INFO/AO)>2 || INFO/DP <= 10' -Oz -o vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.part${SLURM_ARRAY_TASK_ID}.vcf.gz vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.part${SLURM_ARRAY_TASK_ID}.vcf.gz && tabix -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.part${SLURM_ARRAY_TASK_ID}.vcf.gz
bcftools norm -f ${INDEX} vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.part${SLURM_ARRAY_TASK_ID}.vcf.gz | bcftools view -v snps - -Oz -o vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.part${SLURM_ARRAY_TASK_ID}.vcf.gz
tabix -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.part${SLURM_ARRAY_TASK_ID}.vcf.gz
python filter_vcf.py --in vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.part${SLURM_ARRAY_TASK_ID}.vcf.gz --bam reference/read_align/bam/Bf_P.mem2..md.uniq.bam --bam reference/read_align/bam/Bf_M.mem2..md.uniq.bam --out vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.part${SLURM_ARRAY_TASK_ID}.vcf

cat vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.part${SLURM_ARRAY_TASK_ID}.vcf | grep -v "#" |  grep  "TYPE=snp" | \
mawk -F "\t" 'NR==FNR { 
    a[$1] = 1
    next
} 
    $1 in a { 
    print 
}' FNR/callable_genome/tmp/Bf-${id}.contig_inheritance_v0.8.postprocess.list - > result/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.pass.part${SLURM_ARRAY_TASK_ID}.vcf


endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
