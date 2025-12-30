#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=3PZM_calling
#SBATCH --output=log/3PZM_calling_%j_%a.out

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
# count_mutation=${1}
prefix=Bf-${id}
threads=1
INDEX=${user_home}/Amphioxus/new_align/reassembly/read_align/reference/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa
freebayes_folder=${user_home}/Amphioxus/new_align/reassembly/read_align/freebayes_update1
startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`

source ${conda_home}/etc/profile.d/conda.sh
conda activate base

bcftools view  \
    -R ${freebayes_folder}/tmp/${prefix}.out/sites.list -c 1  \
    ${freebayes_folder}/vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz -Ov | \
    bcftools norm -m -any | \
    bcftools norm -a -f ${INDEX} | \
    bcftools view  \
    -v snps -e 'INFO/DP <= 10 || INFO/AO<=5' -Ov  | \
    bgzip -c > vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.het_snp.vcf.gz \
    && tabix -f -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.het_snp.vcf.gz
python script/filter_vcf.py --in vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.het_snp.vcf.gz --bam ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/Bf_P.mem2..md.uniq.bam --bam ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/Bf_M.mem2..md.uniq.bam --out vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.het_snp.filtered.vcf

bcftools view  \
    -R ${freebayes_folder}/tmp/${prefix}.unphased.out/sites.list -c 1 \
    ${freebayes_folder}/vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.unphased.vcf.gz -Ov | \
    bcftools norm -m -any | \
    bcftools norm -a -f ${INDEX} | \
    bcftools view  \
    -v snps -e 'INFO/DP <= 10 || INFO/AO<=5' -Ov | \
    bgzip -c > vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.het_snp.vcf.gz \
    && tabix -f -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.het_snp.vcf.gz
python script/filter_vcf.py --in vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.het_snp.vcf.gz --bam ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/Bf_P.mem2..md.uniq.bam --bam ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/Bf_M.mem2..md.uniq.bam --out vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.het_snp.filtered.vcf

bcftools concat \
    vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.het_snp.filtered.vcf \
    vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.het_snp.filtered.vcf \
    -Ov -o vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.het_snp.filtered.concat_unphased.vcf

python script/extract_alt_ad_minDP2.py --input vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.het_snp.filtered.concat_unphased.vcf --output tmp/${prefix}.alt_ad

awk -F "\t" 'NR==FNR {a[$1,$2]=1;next} {if(!(($1,$2) in a)) {print;}}'  repeatVariants.list tmp/${prefix}.alt_ad.tsv  > tmp/${prefix}.alt_ad.filterrepeat.tsv

python script/filter_alt_allele.py --in tmp/${prefix}.alt_ad.filterrepeat.tsv --bam ${SLURM_ARRAY_TASK_ID} --out result/${prefix}.alt_ad.filterrepeat.filteralt_allsamples.tsv

cat ${user_home}/Amphioxus/new_align/reassembly/read_align/freebayes_update1/filter/indel/indel_regions_${prefix}.bed ${user_home}/Amphioxus/new_align/reassembly/read_align/freebayes_update1/filter/parents_depth/filtered.list | awk -F "\t" 'NR==FNR {a[$1"\t"$3]=1;next;} !($1"\t"$2 in a){OFS="\t";print;}' - result/${prefix}.alt_ad.filterrepeat.filteralt_allsamples.tsv | awk -F "\t" '$7>5' > result/${prefix}.alt_ad.filterrepeat.filteralt_allsamples_removedindel.tsv


endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
