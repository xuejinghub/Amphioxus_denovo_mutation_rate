#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=20G
#SBATCH --job-name=postzygoticDNM
#SBATCH --output=log/postzygoticDNM_%a.out

if [ ! -d tmp ]
then
    mkdir tmp
fi
if [ ! -d vcf ]
then    
    mkdir vcf
fi
if [ ! -d result ]
then
    mkdir result
fi

prefix=Bf-${SLURM_ARRAY_TASK_ID}
threads=2
INDEX=/public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/reference/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa

zcat vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.vcf.gz | grep -v "INDEL" | bcftools view -i 'GT="het"' -Ov - -o vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.het_snp.vcf
cat vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.filtered.het_snp.vcf | grep -v "#" | while read line
do
    contig=$(echo $line | awk -F "\t" '{print $1}')
    pos=$(echo $line | awk -F "\t" '{print $2}')
    ref=$(echo $line | awk -F "\t" '{print $4}')
    filte=$contig":"$pos
    grep "$filter" /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/reference/genome_alignment/minimap2/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.asm20_notremoveSuppl.sort.uniq.bam.contigs0.4.sitelist | awk -F "\t" '{contig="${contig}";if(contig ~ /_P$/) {for(i=4;i<=5;i++){if($i!~/^NA/) {print $i}}} else if(contig ~ /_M$/) {for(i=6;i<=7;i++){if($i!~/^NA/) {print $i}}}}' -  | cut -d \: -f1,2  | while read j; do grep "$j" vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.filtered.het_snp.vcf; done
done

mawk -F "\t" 'NR==FNR { 
    a[$0] = 1
    next
} 
    /^#/ { 
    print
    next
}
    $1 in a { 
    print 
}' inheritance_filter/pass_list/${prefix}.passlist vcf/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.het_snp.vcf > result/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.het_snp.pass.vcf

python filter_vcf.py --in result/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.het_snp.pass.vcf --bam bam/Bf_P.mem2..md.uniq.bam --bam bam/Bf_M.mem2..md.uniq.bam --out result/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.het_snp.pass.bamfiltered.vcf

python extract_alt_ad_minDP2.py   --input result/Bf-${SLURM_ARRAY_TASK_ID}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.bcftools.noMP.het_snp.pass.bamfiltered.vcf   --output updata_result/Bf-${SLURM_ARRAY_TASK_ID}.alt_ad

awk -F "\t" 'NR==FNR {a[$1,$2]=1;next} {if(!(($1,$2) in a)) {print;}}' repeat.list updata_result/${prefix}.alt_ad.tsv  > updata_result/${prefix}.alt_ad.filterrepeat.tsv

cat updata_result/${prefix}.alt_ad.filterrepeat.tsv | cut -f1,2 | while read line; do grep "$line" all_sample_variants.unfiltered.list | awk '$3!="'${prefix}'"' | awk -F "\t" '!array[$1,$2]++' | cut -f1,2; done > tmp/${prefix}_filtered.list

grep -v -f tmp/${prefix}_filtered.list updata_result/${prefix}.alt_ad.filterrepeat.tsv > updata_result/${prefix}.alt_ad.filterrepeat.filteroverlap.tsv

python filter_alt_allele.py --in updata_result/${prefix}.alt_ad.filterrepeat.filteroverlap.tsv --bam ${SLURM_ARRAY_TASK_ID} --out updata_result/${prefix}.alt_ad.filterrepeat.filteroverlap.filteralt_allsamples.tsv