#!/bin/bash
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=freebayes
#SBATCH --output=log/freebayes_%j_%a.out

if [ ! -d filter ]
then
    mkdir filter
fi

cat filter/indel/Bf-*.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.hetero_allele.norm.removedindel.list | awk -F "\t" '!array[$0]++' > filter/indel/Bf-all.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.hetero_allele.norm.removedindel.uniq.list

cat filter/hetero_allele/Bf-*.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.hetero_allele.norm.removedup.vcf |  awk -F "\t" 'FNR==NR {a[$2"\t"$3]=1;next;}  !($2"\t"$3 in a) {OFS="\t";print $0}' - gc/filter_gc/Bf-allsample.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.filtergc.vcf | \
awk -F'\t' '{
    split($2, parts, "_"); 
    for(i in parts) { 
        if(match(parts[i], /^len[0-9]+/)) { 
            len=substr(parts[i], 4)
            break 
        } 
    }
    if($3>100 && (len-$3)>100) {
        print $0
    }
}'  |\
awk -F "\t" 'FNR==NR {a[$4"\t"$1"\t"$2]=1;next;} !($1"\t"$2"\t"$3 in a){OFS="\t";print;}' filter/indel/Bf-all.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.hetero_allele.norm.removedindel.uniq.list - > filter/Bf-allsample.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.filtergc.filterparental_allele_len_indel.vcf

cat filter/Bf-allsample.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.filtergc.filterparental_allele_len_indel.vcf | awk -F "\t" '{print $2"\t"$3-1"\t"$3}' | bedtools intersect -a - -b filter/parents_depth/Bf_P.mem2..md.uniq.bam.depth.per-base.lt5dep.sorted.merged.inherited.bed filter/parents_depth/Bf_M.mem2..md.uniq.bam.depth.per-base.lt5dep.sorted.merged.inherited.bed > filter/parents_depth/filtered.list

awk -F "\t" 'FNR==NR {a[$1"\t"$3]=1;next;} !($2"\t"$3 in a){OFS="\t";print;}' filter/parents_depth/filtered.list filter/Bf-allsample.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.filtergc.filterparental_allele_len_indel.vcf  > filter/Bf-allsample.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.filtergc.filterparentalallele_len_indel_parentsdepth.vcf

awk -F "\t" 'FNR==NR {a[$1"\t"$2"\t"$3]=1;next;} !($1"\t"$2"\t"$3 in a){OFS="\t";print;}' filter/Bf-allsample.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.filtergc.filterparentalallele_len_indel_parentsdepth.vcf final_result/DNM_list.removegc.addsimilar.tsv | wc -l