#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=100G
#SBATCH --job-name=filter_truenonbubble
#SBATCH --output=log/filter_truenonbubble_%j.out

chr=${1}
het=$3
parent=$2
depth=40
species=bf
prefix=${chr}_sim_het${het}_${parent}_dep${depth}_batch${4}
genome_size=$(cat ${user_home}/Amphioxus/bf.dnm/reference/bf.fa.fai | awk -v chr="${chr}" '$1==chr{print $2}')
cpu=20
fastq_dictory=wd/${prefix}/fq
min_score=200
align_software=mem2
parameters=${chr}_sim_het${het}_${parent}_dep${depth}_allPhasedScaffold.min150bp
threads=${cpu}

source /public/home/xuejing/mambaforge/etc/profile.d/conda.sh
conda activate base

mkdir -p wd/${prefix}

cd wd/${prefix}

min_score=200
min_id=0.7
min_cov=0.7
para="min_id"${min_id}_"min_cov"${min_cov}

refrence=fa/${chr}_sim_het${het}_${parent}_dep${depth}_poisson_diploid.fa

if [ ! -d align ]; then
    mkdir align
fi

if [ ! -d result/${para} ]; then
    mkdir -p result/${para}
fi

wd=result/${para}

INDEX=fa/${chr}_sim_het${het}_${parent}_dep${depth}_allPhasedScaffold.min150.fa

cat bam/${prefix}.${align_software}.${parameters}.md.uniq.bam.depth.mosdepth.summary.txt | \
    awk -F "\t" '
    $4>50 && $1!="total" && $4<100 && $1~/^non/ {print $1}
    ' > result/${prefix}.true_nonbubble.list

mawk -F "\t" 'NR==FNR {
    a[$0]=1
    next
}
    ($1 in a) || ($6 in a){
    OFS="\t"
    print $0
    }
' result/${prefix}.true_nonbubble.list align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.paf > align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg3_true_nonbubble.paf

sort -k1,1 -k6,6 align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg3_true_nonbubble.paf > align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg3_true_nonbubble.sort.paf

python3 scripts/paf_filter_rbh.py \
  --paf align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg3_true_nonbubble.sort.paf \
  --min-id $min_id \
  --min-cov $min_cov \
  --out-pairs ${wd}/filtered_pairs_true_nonbubble.${min_id}.${min_cov}.tsv \
  --out-failed-filter ${wd}/failed_pairs_true_nonbubble.${min_id}.${min_cov}.tsv \
  --out-rbh ${wd}/rbh_true_nonbubble.${min_id}.${min_cov}.tsv \
  --rbh-window 100 \
  --end-relax 100 \
  --out-failed-rbh ${wd}/failed_rbh_true_nonbubble.${min_id}.${min_cov}.tsv

cat result/${para}/rbh_true_nonbubble.0.7.0.7.tsv | sed '1d' | awk -F "\t" '{print $1"\n"$5}' | awk -F "\t" 'NR==FNR{a[$1]=1;next;} !($0 in a){print $0}' - result/${prefix}.true_nonbubble.list > result/${para}/rbh_true_nonbubble.0.7.0.7.filtered.tsv

refrence=${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/align_readdepth/sim/pisson_assebly/sim_genome/${prefix}/${chr}_sim_het${het}_${parent}_dep${depth}_hap1.fa

minimap2 -ax asm20 -t 20 -m ${min_score} --secondary=no ${refrence} ${INDEX} | \
    samtools sort -o align/${prefix}to$(basename ${refrence} .fa).sim.sort.bam - 
samtools index align/${prefix}to$(basename ${refrence} .fa).sim.sort.bam

minimap2 -cx asm20 -t 20 -m ${min_score} --secondary=no ${refrence} ${INDEX} > align/${prefix}to$(basename ${refrence} .fa).sim.paf

realnon_len=$(cat result/${para}/rbh_true_nonbubble.0.7.0.7.filtered.tsv | awk -F '_' '{for(i=1;i<=NF;i++)if($i~/^len/){split($i,a,"n");sum+=a[2];break}}END{print sum}')

cat sim_genome/${prefix}/nonvariant_regions.bed | awk -F "\t" '{OFS="\t";print $1,$3-1,$3,"region_"NR}' | bedtools closest -t all -a - -b sim_genome/${prefix}/genome_simulated_variants.vcf.gz -k 1 -io -iu -D ref > sim_genome/${prefix}/nonvariant_regions_end.bed

cat sim_genome/${prefix}/nonvariant_regions.bed | awk -F "\t" '{OFS="\t";print $1,$3-1,$3,"region_"NR}' | bedtools closest -t all -a - -b sim_genome/${prefix}/genome_simulated_variants.vcf.gz -k 1 -io -id -D ref  > sim_genome/${prefix}/nonvariant_regions_start.bed

awk -F "\t" 'NR==FNR{a[$4]=$6;next;} $4 in a{OFS="\t";print $1,a[$4],$6,$4}'  sim_genome/${prefix}/nonvariant_regions_start.bed sim_genome/${prefix}/nonvariant_regions_end.bed > sim_genome/${prefix}/nonvariant_regions_expand2nonvariant.bed

cat align/${prefix}to$(basename ${refrence} .fa).sim.paf  | sed 's/_hap1//g' |  awk -F "\t" 'NR==FNR{a[$1]=1;next;} $1 in a{OFS="\t";print $6,$8,$9,$1,$2,$10/$2}' result/${para}/rbh_true_nonbubble.0.7.0.7.filtered.tsv  - | bedtools coverage -hist -a stdin -b ${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/align_readdepth/sim/pisson_assebly/sim_genome/${prefix}/nonvariant_regions_expand2nonvariant.bed | tail | awk -F "\t" -v realnon_len="${realnon_len}"  -v para="${para}" 'BEGIN{OFS="\t";print "para","align_base","nonanlign_base","all_realnon_base","False_positive"}{if($1=="all" && $2==0){nonalign_base=$3}}END{OFS="\t";print para,realnon_len-nonalign_base,nonalign_base,realnon_len,nonalign_base/realnon_len}' > result/${para}/rbh_true_nonbubble.0.7.0.7.checked.stats

calculate_total_length() {
    local tsv_file=$1
    t_col=$(head -1 $tsv_file | awk -F "\t" '{for(i=1;i<=NF;i++)if($i="t"){print i;break}}')
    cat "$tsv_file" | sed '1d' | cut -f1,${t_col} | awk '{print $1;print $2}' | sort | uniq | awk -F '_' '{for(i=1;i<=NF;i++)if($i~/^len/){split($i,a,"n");sum+=a[2];break}}END{print sum}'
}
calculate_total_length_truenonbubble() {
    local tsv_file=$1
    t_col=$(head -1 $tsv_file | awk -F "\t" '{for(i=1;i<=NF;i++)if($i="t"){print i;break}}')
    cat "$tsv_file" | sed '1d' | cut -f1,${t_col} | awk '{print $1;print $2}' | sort | uniq  | awk -F "\t" 'NR==FNR {a[$1]=1;next} !($1 in a){print $1}' - result/${prefix}.true_nonbubble.list  | awk -F '_' '{for(i=1;i<=NF;i++)if($i~/^len/){split($i,a,"n");sum+=a[2];break}}END{print sum}'
}

calculate_contig_num() {
    local tsv_file=$1
    t_col=$(head -1 $tsv_file | awk -F "\t" '{for(i=1;i<=NF;i++)if($i="t"){print i;break}}')
    cat "$tsv_file" | sed '1d' | cut -f1,${t_col} | awk '{print $1;print $2}' | sort | uniq | wc -l
}

calculate_contig_num_truenonbubble() {
    local tsv_file=$1
    t_col=$(head -1 $tsv_file | awk -F "\t" '{for(i=1;i<=NF;i++)if($i="t"){print i;break}}')
    cat "$tsv_file" | sed '1d' | cut -f1,${t_col} | awk '{print $1;print $2}' | sort | uniq | awk -F "\t" 'NR==FNR {a[$1]=1;next} !($1 in a){print $1}'  - result/${prefix}.true_nonbubble.list | wc -l
}

calculate_contig_num_truenonbubble_list() {
    local tsv_file=$1
    t_col=$(head -1 $tsv_file | awk -F "\t" '{for(i=1;i<=NF;i++)if($i="t"){print i;break}}')
    cat "$tsv_file" | sed '1d' | cut -f1,${t_col} | awk '{print $1;print $2}' | sort | uniq | awk -F "\t" 'NR==FNR {a[$1]=1;next} !($1 in a){print $1}' - result/${prefix}.true_nonbubble.list > ${prefix}.true_nonbubble.checked.list
}

total_length_bubblecontigs=$(calculate_total_length ${wd}/rbh_bubblecontigs.${min_id}.${min_cov}.tsv)
total_length_bubblelikecontigs=$(calculate_total_length ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.tsv)
total_length_true_nonbubble=$(calculate_total_length_truenonbubble ${wd}/rbh_true_nonbubble.${min_id}.${min_cov}.tsv)
echo "Total length of real nonbubble contigs: $total_length_true_nonbubble"
echo "Total length of bubble contigs: $total_length_bubblecontigs"
echo "Total length of bubble like contigs: $total_length_bubblelikecontigs"

total_length_filtered_bubblecontigs=$(calculate_total_length ${wd}/filtered_pairs_bubblecontigs.${min_id}.${min_cov}.tsv)
total_length_filtered_bubblelikecontigs=$(calculate_total_length ${wd}/filtered_pairs_bubblelikecontigs.${min_id}.${min_cov}.tsv)
total_length_filtered_true_nonbubble=$(calculate_total_length_truenonbubble ${wd}/filtered_pairs_true_nonbubble.${min_id}.${min_cov}.tsv)
echo "Total length of filtered real nonbubble contigs: $total_length_filtered_true_nonbubble"
echo "Total length of filtered bubble contigs: $total_length_filtered_bubblecontigs"
echo "Total length of filtered bubble like contigs: $total_length_filtered_bubblelikecontigs"

total_length_rbh_failed_bubblecontigs=$(calculate_total_length ${wd}/failed_rbh_bubblecontigs.${min_id}.${min_cov}.tsv)
total_length_rbh_failed_bubblelikecontigs=$(calculate_total_length ${wd}/failed_rbh_bubblelikecontigs.${min_id}.${min_cov}.tsv)
total_length_rbh_failed_true_nonbubble=$(calculate_total_length_truenonbubble ${wd}/failed_rbh_true_nonbubble.${min_id}.${min_cov}.tsv)
echo "Total length of failed rbh real nonbubble contigs: $total_length_rbh_failed_true_nonbubble"
echo "Total length of failed rbh bubble contigs: $total_length_rbh_failed_bubblecontigs"
echo "Total length of failed rbh bubble like contigs: $total_length_rbh_failed_bubblelikecontigs"


contig_num_bubblecontigs=$(calculate_contig_num ${wd}/rbh_bubblecontigs.${min_id}.${min_cov}.tsv)
contig_num_bubblelikecontigs=$(calculate_contig_num ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.tsv)
contig_num_true_nonbubble=$(calculate_contig_num_truenonbubble ${wd}/rbh_true_nonbubble.${min_id}.${min_cov}.tsv)
echo "Total number of real nonbubble contigs: $contig_num_true_nonbubble"
echo "Total number of bubble contigs: $contig_num_bubblecontigs"
echo "Total number of bubble like contigs: $contig_num_bubblelikecontigs"

filtered_contig_num_bubblecontigs=$(calculate_contig_num ${wd}/filtered_pairs_bubblecontigs.${min_id}.${min_cov}.tsv)
filtered_contig_num_bubblelikecontigs=$(calculate_contig_num ${wd}/filtered_pairs_bubblelikecontigs.${min_id}.${min_cov}.tsv)
filtered_contig_num_true_nonbubble=$(calculate_contig_num_truenonbubble ${wd}/filtered_pairs_true_nonbubble.${min_id}.${min_cov}.tsv)
echo "Total number of filtered real nonbubble contigs: $filtered_contig_num_true_nonbubble"
echo "Total number of filtered bubble contigs: $filtered_contig_num_bubblecontigs"
echo "Total number of filtered bubble like contigs: $filtered_contig_num_bubblelikecontigs"

rbh_failed_contig_num_bubblecontigs=$(calculate_contig_num ${wd}/failed_rbh_bubblecontigs.${min_id}.${min_cov}.tsv)
rbh_failed_contig_num_bubblelikecontigs=$(calculate_contig_num ${wd}/failed_rbh_bubblelikecontigs.${min_id}.${min_cov}.tsv)
rbh_failed_contig_num_true_nonbubble=$(calculate_contig_num_truenonbubble ${wd}/failed_rbh_true_nonbubble.${min_id}.${min_cov}.tsv)
echo "Total number of failed rbh real nonbubble contigs: $rbh_failed_contig_num_true_nonbubble"
echo "Total number of failed rbh bubble contigs: $rbh_failed_contig_num_bubblecontigs"
echo "Total number of failed rbh bubble like contigs: $rbh_failed_contig_num_bubblelikecontigs"


## Output statistic std
awk -F "\t" \
    -v para=$para \
    -v contig_num_bubblecontigs=$contig_num_bubblecontigs \
    -v contig_num_bubblelikecontigs=$contig_num_bubblelikecontigs \
    -v filtered_contig_num_bubblecontigs=$filtered_contig_num_bubblecontigs \
    -v filtered_contig_num_bubblelikecontigs=$filtered_contig_num_bubblelikecontigs \
    -v rbh_failed_contig_num_bubblecontigs=$rbh_failed_contig_num_bubblecontigs \
    -v rbh_failed_contig_num_bubblelikecontigs=$rbh_failed_contig_num_bubblelikecontigs \
    -v total_length_bubblecontigs=$total_length_bubblecontigs \
    -v total_length_bubblelikecontigs=$total_length_bubblelikecontigs \
    -v total_length_filtered_bubblecontigs=$total_length_filtered_bubblecontigs \
    -v total_length_filtered_bubblelikecontigs=$total_length_filtered_bubblelikecontigs \
    -v total_length_rbh_failed_bubblecontigs=$total_length_rbh_failed_bubblecontigs \
    -v total_length_rbh_failed_bubblelikecontigs=$total_length_rbh_failed_bubblelikecontigs \
    -v total_length_true_nonbubble=$total_length_true_nonbubble \
    -v total_length_filtered_true_nonbubble=$total_length_filtered_true_nonbubble \
    -v total_length_rbh_failed_true_nonbubble=$total_length_rbh_failed_true_nonbubble \
    -v contig_num_true_nonbubble=$contig_num_true_nonbubble \
    -v filtered_contig_num_true_nonbubble=$filtered_contig_num_true_nonbubble \
    -v rbh_failed_contig_num_true_nonbubble=$rbh_failed_contig_num_true_nonbubble \
    'BEGIN{
        print         "para\t#final_bubblecontigs\t#filtered_bubblecontigs\t#rbh_failed_bubblecontigs\t#final_bubblelikecontigs\t#filtered_bubblelikecontigs\t#rbh_failed_bubblelikecontigs\t#final_true_nonbubblecontigs\t#filtered_true_nonbubblecontigs\t#rbh_failed_true_nonbubblecontigs\tlen_final_bubblecontigs\tlen_filtered_bubblecontigs\tlen_rbh_failed_bubblecontigs\tlen_final_bubblelikecontigs\tlen_filtered_bubblelikecontigs\tlen_rbh_failed_bubblelikecontigs\tlen_final_true_nonbubblecontigs\tlen_filtered_true_nonbubblecontigs\tlen_rbh_failed_true_nonbubblecontigs"
    }END{ 
        print para"\t"contig_num_bubblecontigs"\t"filtered_contig_num_bubblecontigs"\t"rbh_failed_contig_num_bubblecontigs"\t"contig_num_bubblelikecontigs"\t"filtered_contig_num_bubblelikecontigs"\t"rbh_failed_contig_num_bubblelikecontigs"\t"contig_num_true_nonbubble"\t"filtered_contig_num_true_nonbubble"\t"rbh_failed_contig_num_true_nonbubble"\t"total_length_bubblecontigs"\t"total_length_filtered_bubblecontigs"\t"total_length_rbh_failed_bubblecontigs"\t"total_length_bubblelikecontigs"\t"total_length_filtered_bubblelikecontigs"\t"total_length_rbh_failed_bubblelikecontigs"\t"total_length_true_nonbubble"\t"total_length_filtered_true_nonbubble"\t"total_length_rbh_failed_true_nonbubble
    }' > ${wd}/${para}.statistic2.tsv

calculate_contig_num_truenonbubble_list ${wd}/rbh_true_nonbubble.${min_id}.${min_cov}.tsv 



