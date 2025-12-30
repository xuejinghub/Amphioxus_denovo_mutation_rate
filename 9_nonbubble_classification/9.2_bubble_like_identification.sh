#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=100G
#SBATCH --job-name=align_readdepth
#SBATCH --output=log/align_readdepth_%j.out

prefix=Bf_${1}

mkdir -p wd/${prefix}

cd wd/${prefix}

min_score=200
min_id=$2
min_cov=$3
para="min_id"${min_id}_"min_cov"${min_cov}

if [ ! -d align ]; then
    mkdir align
fi

if [ ! -d result/${para} ]; then
    mkdir -p result/${para}
fi

wd=result/${para}

INDEX=${user_home}/Amphioxus/new_align/reassembly/read_align/reference/${prefix}_platanus_i3_allPhasedScaffold.rename.min150bp.fa

sort -k1,1 -k6,6 align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg2_bubblecontigs.paf > align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg2_bubblecontigs.sort.paf

python3 script/paf_filter_rbh.py \
  --paf align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg2_bubblecontigs.sort.paf \
  --min-id $min_id \
  --min-cov $min_cov \
  --out-pairs ${wd}/filtered_pairs_bubblecontigs.${min_id}.${min_cov}.tsv \
  --out-failed-filter ${wd}/failed_pairs_bubblecontigs.${min_id}.${min_cov}.tsv \
  --out-rbh ${wd}/rbh_bubblecontigs.${min_id}.${min_cov}.tsv \
  --rbh-window 100 \
  --end-relax 100 \
  --out-failed-rbh ${wd}/failed_rbh_bubblecontigs.${min_id}.${min_cov}.tsv

sort -k1,1 -k6,6 align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg1.paf > align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg1.sort.paf  
python3 script/paf_filter_rbh.py \
  --paf align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg1.sort.paf \
  --min-id $min_id \
  --min-cov $min_cov \
  --out-pairs ${wd}/filtered_pairs_bubblelikecontigs.${min_id}.${min_cov}.tsv \
  --out-failed-filter ${wd}/failed_pairs_bubblelikecontigs.${min_id}.${min_cov}.tsv \
  --out-rbh ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.tsv \
  --rbh-window 100 \
  --end-relax 100 \
  --out-failed-rbh ${wd}/failed_rbh_bubblelikecontigs.${min_id}.${min_cov}.tsv

calculate_total_length() {
    local tsv_file=$1
    t_col=$(head -1 $tsv_file | awk -F "\t" '{for(i=1;i<=NF;i++)if($i="t"){print i;break}}')
    cat "$tsv_file" | sed '1d' | cut -f1,${t_col} | awk '{print $1;print $2}' | sort | uniq | awk -F '_' '{for(i=1;i<=NF;i++)if($i~/^len/){split($i,a,"n");sum+=a[2];break}}END{print sum}'
}

calculate_contig_num() {
    local tsv_file=$1
    t_col=$(head -1 $tsv_file | awk -F "\t" '{for(i=1;i<=NF;i++)if($i="t"){print i;break}}')
    cat "$tsv_file" | sed '1d' | cut -f1,${t_col} | awk '{print $1;print $2}' | sort | uniq | wc -l
}

## statistics of output
total_length_bubblecontigs=$(calculate_total_length ${wd}/rbh_bubblecontigs.${min_id}.${min_cov}.tsv)
total_length_bubblelikecontigs=$(calculate_total_length ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.tsv)
echo "Total length of bubble contigs: $total_length_bubblecontigs"
echo "Total length of bubble like contigs: $total_length_bubblelikecontigs"

total_length_filtered_bubblecontigs=$(calculate_total_length ${wd}/filtered_pairs_bubblecontigs.${min_id}.${min_cov}.tsv)
total_length_filtered_bubblelikecontigs=$(calculate_total_length ${wd}/filtered_pairs_bubblelikecontigs.${min_id}.${min_cov}.tsv)
echo "Total length of filtered bubble contigs: $total_length_filtered_bubblecontigs"
echo "Total length of filtered bubble like contigs: $total_length_filtered_bubblelikecontigs"

total_length_rbh_failed_bubblecontigs=$(calculate_total_length ${wd}/failed_rbh_bubblecontigs.${min_id}.${min_cov}.tsv)
total_length_rbh_failed_bubblelikecontigs=$(calculate_total_length ${wd}/failed_rbh_bubblelikecontigs.${min_id}.${min_cov}.tsv)
echo "Total length of failed rbh bubble contigs: $total_length_rbh_failed_bubblecontigs"
echo "Total length of failed rbh bubble like contigs: $total_length_rbh_failed_bubblelikecontigs"

contig_num_bubblecontigs=$(calculate_contig_num ${wd}/rbh_bubblecontigs.${min_id}.${min_cov}.tsv)
contig_num_bubblelikecontigs=$(calculate_contig_num ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.tsv)
echo "Total number of bubble contigs: $contig_num_bubblecontigs"
echo "Total number of bubble like contigs: $contig_num_bubblelikecontigs"

filtered_contig_num_bubblecontigs=$(calculate_contig_num ${wd}/filtered_pairs_bubblecontigs.${min_id}.${min_cov}.tsv)
filtered_contig_num_bubblelikecontigs=$(calculate_contig_num ${wd}/filtered_pairs_bubblelikecontigs.${min_id}.${min_cov}.tsv)
echo "Total number of filtered bubble contigs: $filtered_contig_num_bubblecontigs"
echo "Total number of filtered bubble like contigs: $filtered_contig_num_bubblelikecontigs"

rbh_failed_contig_num_bubblecontigs=$(calculate_contig_num ${wd}/failed_rbh_bubblecontigs.${min_id}.${min_cov}.tsv)
rbh_failed_contig_num_bubblelikecontigs=$(calculate_contig_num ${wd}/failed_rbh_bubblelikecontigs.${min_id}.${min_cov}.tsv)
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
    'BEGIN{
        print         "para\t#final_bubblecontigs\t#filtered_bubblecontigs\t#rbh_failed_bubblecontigs\t#final_bubblelikecontigs\t#filtered_bubblelikecontigs\t#rbh_failed_bubblelikecontigs\tlen_final_bubblecontigs\tlen_filtered_bubblecontigs\tlen_rbh_failed_bubblecontigs\tlen_final_bubblelikecontigs\tlen_filtered_bubblelikecontigs\tlen_rbh_failed_bubblelikecontigs"
    }END{ 
        print para"\t"contig_num_bubblecontigs"\t"filtered_contig_num_bubblecontigs"\t"rbh_failed_contig_num_bubblecontigs"\t"contig_num_bubblelikecontigs"\t"filtered_contig_num_bubblelikecontigs"\t"rbh_failed_contig_num_bubblelikecontigs"\t"total_length_bubblecontigs"\t"total_length_filtered_bubblecontigs"\t"total_length_rbh_failed_bubblecontigs"\t"total_length_bubblelikecontigs"\t"total_length_filtered_bubblelikecontigs"\t"total_length_rbh_failed_bubblelikecontigs
    }' > ${wd}/${para}.statistic.tsv


