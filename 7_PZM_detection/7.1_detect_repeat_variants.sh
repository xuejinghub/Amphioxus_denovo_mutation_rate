#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=repeatVariants
#SBATCH --output=log/repeatVariants_%j_%a.out

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
INDEX=${user_home}/Amphioxus/new_align/reassembly/read_align/reference/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa
freebayes_folder=${user_home}/Amphioxus/new_align/reassembly/read_align/freebayes_update1/

bcftools concat ${freebayes_folder}/vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.unphased.vcf.gz \
    ${freebayes_folder}/vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz | \
    bcftools norm -m -any | \
    bcftools norm -a -f ${INDEX} -Ov | \
    awk -F "\t" '$0!~/^#/ {print $1"\t"$2}' > tmp/${prefix}.pos

# cat tmp/Bf-*.pos | awk -F "\t" 'array[$0]++' | awk -F "\t" '!array[$0]++' > repeatVariants.list