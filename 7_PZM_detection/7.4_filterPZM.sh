#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=freebayes
#SBATCH --output=log/freebayes_%j_%a.out

if [ ! -d filter/parents_depth ]
then
    mkdir filter/parents_depth
fi
cat result/Bf-*.alt_ad.filterrepeat.filteralt_allsamples_removedindel.tsv | grep -v "CHROM" | cat <(head -1 result/Bf-106.alt_ad.filterrepeat.filteralt_allsamples_removedindel.tsv) - > result/Bf-allsample.alt_ad.filterrepeat.filteralt_allsamples_removedindel1.tsv

cat result/Bf-allsample.alt_ad.filterrepeat.filteralt_allsamples_removedindel1.tsv | sed '1d' | awk -F "\t" '{print $1"\t"$2-1"\t"$2}' | bedtools intersect -a - -b /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/freebayes/filter/parents_depth/Bf_P.mem2..md.uniq.bam.depth.per-base.lt5dep.sorted.merged.inherited.bed /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/freebayes/filter/parents_depth/Bf_M.mem2..md.uniq.bam.depth.per-base.lt5dep.sorted.merged.inherited.bed > filter/parents_depth/filtered.list

awk -F "\t" 'FNR==NR {a[$1"\t"$3]=1;next;} !($1"\t"$2 in a){OFS="\t";print;}' filter/parents_depth/filtered.list result/Bf-allsample.alt_ad.filterrepeat.filteralt_allsamples_removedindel1.tsv > result/Bf-allsample.alt_ad.filterrepeat.filteralt_allsamples_removedindel_parentsdepth.tsv

