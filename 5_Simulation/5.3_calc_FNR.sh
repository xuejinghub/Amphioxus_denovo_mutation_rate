#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=stats
#SBATCH --output=log/stats_%j_%a.out

stats_function() {
    local id=$1
    local mut=$2
    bin_count=$(ls FNR/result/Bf${id}_mut${mut}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.pass.part*.vcf | wc -l)
    cat FNR/result/Bf${id}_mut${mut}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.pass.part*.vcf > Bf${id}_mut${mut}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.pass.vcf
    FN=$(awk -F "\t" 'FNR==NR {a[$1"\t"$2]=1;next} !($1"\t"$2 in a) {OFS="\t";print;}' FNR/Bf${id}_mut${mut}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.pass.vcf simulate_v3.0/art_v3.0FNR/Bf${id}_mut${mut}.mut.txt | awk -F "\t" 'FNR==NR {a[$1]=1;next} $1 in a {OFS="\t";print;}' FNR/callable_genome/tmp/Bf-${id}.contig_inheritance_v0.8.postprocess.list  - | wc -l)
    FP=$(awk -F "\t" 'FNR==NR {a[$1"\t"$2]=1;next} ($1"\t"$2 in a) {OFS="\t";print;}' FNR/Bf${id}_mut${mut}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.pass.vcf simulate_v3.0/art_v3.0FNR/Bf${id}_mut${mut}.mut.txt | wc -l)
    awk -F "\t" 'FNR==NR {a[$1"\t"$2]=1;next} !($1"\t"$2 in a) {OFS="\t";print;}' FNR/Bf${id}_mut${mut}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.pass.vcf simulate_v3.0/art_v3.0FNR/Bf${id}_mut${mut}.mut.txt | awk -F "\t" 'FNR==NR {a[$1]=1;next} $1 in a {OFS="\t";print;}' FNR/callable_genome/tmp/Bf-${id}.contig_inheritance_v0.8.postprocess.list - > Bf${id}_mut${mut}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.pass.FN.list
    all=$((${FP}+${FN}))
    FNR=$(awk -v FN=${FN} -v all=${all} 'BEGIN{print FN/all*100}')
    echo -e "Bf-"${id}"-mut${mut}\t"$FN"\t"$FP"\t"$bin_count"\t"${FNR}
}

for i in {1..3}
do
    stats_function ${i} 500
    stats_function ${i} 1000
done

