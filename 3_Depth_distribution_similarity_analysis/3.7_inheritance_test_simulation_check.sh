#!/bin/bash
#SBATCH -w node80
#SBATCH -c 10
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=3inh_check
#SBATCH --output=log/3inh_check_%j_%a.out

chr=Chr3
het=0.03
depth=40
species=bf
batch=${1}

if [ ${batch} == 1 ]; then
    other=2
else
    parent=1
fi

prefix=${chr}_sim_het${het}_F1batch${batch}_dep${depth}
genome_size=$(cat ${user_home}/Amphioxus/bf.dnm/reference/bf.fa.fai | awk -v chr="${chr}" '$1==chr{print $2}')
refrence=${user_home}/Amphioxus/bf.dnm/reference/bf.fa
cpu=20
fastq_dictory=${user_home}/Amphioxus/new_align/reassembly/read_align/cosin_similarity/simulation/wd/${prefix}/fq
min_score=200
align_software=mem2
parameters=${chr}_sim_het${het}_MP_dep${depth}_allPhasedScaffold.min150.rename
threads=${cpu}
script_dir=${user_home}/Amphioxus/new_align/reassembly/read_align/cosin_similarity/simulation/script

sim_genome_dir=${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/align_readdepth/sim/pisson_assebly/sim_genome

if [ ! -d wd/${prefix} ]; then
    mkdir -p wd/${prefix}
fi

cd wd/${prefix}

if [ ! -d tmp ]; then
    mkdir tmp
fi

if [ ! -d fq ]; then
    mkdir fq
fi

if [ ! -d fa ]; then
    mkdir fa
fi

if [ ! -d bam ]; then
    mkdir bam
fi

if [ ! -d align ]; then
    mkdir align
fi

if [ ! -d result ]; then
    mkdir result
fi


cat ${sim_genome_dir}/Chr3_sim_het0.03_M_dep40_batch2/${chr}_sim_het${het}_M_dep${depth}_hap1.fa fa/${chr}_sim_het${het}_M_dep${depth}_hap2.fa > fa/${chr}_sim_het${het}_M_dep${depth}_2haps.fa
cat ${sim_genome_dir}/Chr3_sim_het0.03_P_dep40_batch1/${chr}_sim_het${het}_P_dep${depth}_hap1.fa fa/${chr}_sim_het${het}_P_dep${depth}_hap2.fa > fa/${chr}_sim_het${het}_P_dep${depth}_2haps.fa


minimap2 -cx asm20 -t 5 --cs fa/${chr}_sim_het${het}_M_dep${depth}_2haps.fa fa/${chr}_sim_het${het}_M_dep40_allPhasedScaffold.min150.rename.fa > align/${chr}_dep${depth}_contigs_to_chr_sim_het${het}_M.paf
minimap2 -cx asm20 -t 5 --cs fa/${chr}_sim_het${het}_P_dep${depth}_2haps.fa fa/${chr}_sim_het${het}_P_dep40_allPhasedScaffold.min150.rename.fa > align/${chr}_dep${depth}_contigs_to_chr_sim_het${het}_P.paf

cat align/${chr}_dep${depth}_contigs_to_chr_sim_het${het}_P.paf | awk -F "\t" '{if($1 in a){if($12 > a[$1]) {b[$1]=$6} else if($12 == a[$1] && $12==0){b[$1]="all"}} else {a[$1]=$12;b[$1]=$6}}END{for(i in b){print i"\t"b[i]}}' | awk -F "\t" '$2~/hap2/ || $2~/all/' > align/${chr}_dep${depth}_contigs_to_chr_sim_het${het}_P_TruePos.list

cat align/${chr}_dep${depth}_contigs_to_chr_sim_het${het}_M.paf | awk -F "\t" '{if($1 in a){if($12 > a[$1]) {b[$1]=$6} else if($12 == a[$1] && $12==0){b[$1]="all"}} else {a[$1]=$12;b[$1]=$6}}END{for(i in b){print i"\t"b[i]}}' | awk -F "\t" '$2~/hap2/ || $2~/all/' > align/${chr}_dep${depth}_contigs_to_chr_sim_het${het}_M_TruePos.list

cat align/${chr}_dep${depth}_contigs_to_chr_sim_het${het}_P_TruePos.list align/${chr}_dep${depth}_contigs_to_chr_sim_het${het}_M_TruePos.list > align/${chr}_dep${depth}_contigs_to_chr_sim_het${het}_MP_TruePos.list

truePos=$(cat postprocess_result/${chr}_sim_het${het}_F1batch${batch}_dep${depth}.contig_inheritance_v0.8.postprocess_v0.9.out | grep "PASS" | awk -F "\t" 'NR==FNR{a[$1]=1;next;} $1 in a{sum+=1;sum_len+=$17}END{print sum"\t"sum_len}' align/${chr}_dep${depth}_contigs_to_chr_sim_het${het}_MP_TruePos.list -)
Pos=$(cat postprocess_result/${chr}_sim_het${het}_F1batch${batch}_dep${depth}.contig_inheritance_v0.8.postprocess_v0.9.out | grep "PASS" | awk -F "\t" '{sum+=1;sum_len+=$17}END{print sum"\t"sum_len}')
TP=$(cat align/${chr}_dep${depth}_contigs_to_chr_sim_het${het}_MP_TruePos.list |awk -F '_' '{for(i=1;i<=NF;i++)if($i~/^len/){split($i,a,"n");sum+=a[2];break}}END{print NR"\t"sum}')

FP=$(cat postprocess_result/${chr}_sim_het${het}_F1batch${batch}_dep${depth}.contig_inheritance_v0.8.postprocess_v0.9.out | grep "PASS" | awk -F "\t" 'NR==FNR{a[$1]=1;next;} !($1 in a){sum+=1;sum_len+=$17}END{print sum"\t"sum_len}' align/${chr}_dep${depth}_contigs_to_chr_sim_het${het}_MP_TruePos.list -)

echo -e "TruePos\t${truePos}"
echo -e "Pos\t${Pos}"
echo -e "TP\t${TP}"
echo -e "FP\t${FP}"
