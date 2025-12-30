#!/bin/bash
#SBATCH -w node80
#SBATCH -c 2
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=freebayes
#SBATCH --output=log/freebayes_%j_%a.out


prefix=Bf-${SLURM_ARRAY_TASK_ID}_mut1000
fastq_dictory=FNR/simulation_data
align_software=mem2
parameters="Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp"
threads=50
INDEX=reference/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa
startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
## Step 1: Freebayes calling variants
if [ ! -e vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.sim_freebayes.vcf.gz.tbi ]
then
freebayes -f ${INDEX} \
    -m 0 \
    --legacy-gls \
    -t ${CG_dir}/Bf-${SLURM_ARRAY_TASK_ID}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.callable_filtered.remove3foldsites.remove_unknow_phased_contig.bed \
    bam/${prefix}.${align_software}.${parameters}.md.uniq.bam \
    | bgzip -c > vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.sim_freebayes.vcf.gz
# bcftools sort \
#     -Oz \
#     -o vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.sim_freebayes.sorted.vcf.gz \
#     vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.sim_freebayes.vcf.gz
tabix -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.sim_freebayes.vcf.gz
fi

source ${conda_home}/etc/profile.d/conda.sh
conda activate base

## Step 2: Filter variants
bcftools isec \
    -n~100 \
    -c all \
    -p tmp/${prefix}.out \
    vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.sim_freebayes.vcf.gz  \
    ${user_home}/Amphioxus/new_align/reassembly/read_align/freebayes_update1/vcf/Bf_M.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz \
    ${user_home}/Amphioxus/new_align/reassembly/read_align/freebayes_update1/vcf/Bf_P.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz

cut -f1,2 tmp/${prefix}.out/sites.txt > tmp/${prefix}.out/sites.list
bcftools view  \
    -R tmp/${prefix}.out/sites.list \
    vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.sim_freebayes.vcf.gz -Ov | \
    bcftools norm -m -any | \
    bcftools norm -a -f ${INDEX} | \
    bcftools view  \
    -v snps -e '(INFO/DP)-(INFO/AO)>2 || INFO/DP <= 5' \
    -Oz -o vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.vcf.gz - \
    && tabix -f -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.vcf.gz
python ${user_home}/Amphioxus/new_align/reassembly/read_align/freebayes_update1/script/filter_vcf.py --in vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.vcf.gz --bam ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/Bf_P.mem2..md.uniq.bam --bam ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/Bf_M.mem2..md.uniq.bam --out vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.filtered.vcf

cat vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.filtered.vcf | grep -v "#" | awk  '$0!~/TYPE=indel/' > result/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.norm.filtered.vcf

source ${conda_home}/etc/profile.d/conda.sh
conda activate freebayes

# Step 3: Freebayes calling unphased variants 
if [ ! -e vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.sim_freebayes.unphased.vcf.gz.tbi ]
then
freebayes -f ${INDEX} \
    -m 0 \
    --legacy-gls \
    -t ${CG_dir}/Bf-${SLURM_ARRAY_TASK_ID}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.uncallableE3.bed \
    bam/${prefix}.${align_software}.${parameters}.md.uniq.bam \
    | bgzip -c > vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.sim_freebayes.unphased.vcf.gz
tabix -f -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.sim_freebayes.unphased.vcf.gz
fi


source ${conda_home}/etc/profile.d/conda.sh
conda activate base

## Step 4: Filter unphased variants
bcftools isec \
    -n~100 \
    -c all \
    -p tmp/${prefix}.unphased.out \
    vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.sim_freebayes.unphased.vcf.gz  \
    ${user_home}/Amphioxus/new_align/reassembly/read_align/freebayes_update1/vcf/Bf_M.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz \
    ${user_home}/Amphioxus/new_align/reassembly/read_align/freebayes_update1/vcf/Bf_P.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz

cut -f1,2 tmp/${prefix}.unphased.out/sites.txt > tmp/${prefix}.unphased.out/sites.list
bcftools view  \
    -R tmp/${prefix}.unphased.out/sites.list \
    vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.sim_freebayes.unphased.vcf.gz -Ov | \
    bcftools norm -m -any | \
    bcftools norm -a -f ${INDEX} | \
    bcftools view  \
    -v snps -e '(INFO/DP)-(INFO/AO)>2 || INFO/DP <= 5' \
    -Oz -o vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.vcf.gz - \
    && tabix -f -p vcf vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.vcf.gz
python ${user_home}/Amphioxus/new_align/reassembly/read_align/freebayes_update1/script/filter_vcf.py --in vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.vcf.gz --bam ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/Bf_P.mem2..md.uniq.bam --bam ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/Bf_M.mem2..md.uniq.bam --out vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.filtered.vcf

cat vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.filtered.vcf | grep -v "#" |  awk  '$0!~/TYPE=indel/' > result/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.unphased.noMP.norm.filtered.vcf

endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
