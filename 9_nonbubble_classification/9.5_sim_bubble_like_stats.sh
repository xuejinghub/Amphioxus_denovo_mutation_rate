#!/bin/bash
#SBATCH -w node80
#SBATCH -c 20
#SBATCH --time=48:00:00                     
#SBATCH --mem=100G
#SBATCH --job-name=filter_bubble
#SBATCH --output=log/filter_bubble_%j.out

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

source ${conda_home}/etc/profile.d/conda.sh
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
    awk -F "\t" 'BEGIN{OFS="\t";print "interval","total","total_len","#bubble","bubble_len","#nonbubble","nonbubble_len"}
    $4>0 && $1!="total" && $4<=10 {sum10+=1;sum_len10+=$2;if($1!~/^non/){sum_bubble10+=1;sum_bubble_len10+=$2} else {sum_nonbubble10+=1;sum_nonbubble_len10+=$2}} 
    $4>10 && $1!="total" && $4<=20 {sum20+=1;sum_len20+=$2;if($1!~/^non/){sum_bubble20+=1;sum_bubble_len20+=$2} else {sum_nonbubble20+=1;sum_nonbubble_len20+=$2}}
    $4>20 && $1!="total" && $4<=30 {sum30+=1;sum_len30+=$2;if($1!~/^non/){sum_bubble30+=1;sum_bubble_len30+=$2} else {sum_nonbubble30+=1;sum_nonbubble_len30+=$2}}
    $4>30 && $1!="total" && $4<=40 {sum40+=1;sum_len40+=$2;if($1!~/^non/){sum_bubble40+=1;sum_bubble_len40+=$2} else {sum_nonbubble40+=1;sum_nonbubble_len40+=$2}}
    $4>40 && $1!="total" && $4<=50 {sum50+=1;sum_len50+=$2;if($1!~/^non/){sum_bubble50+=1;sum_bubble_len50+=$2} else {sum_nonbubble50+=1;sum_nonbubble_len50+=$2}}
    $4>50 && $1!="total" && $4<=60 {sum60+=1;sum_len60+=$2;if($1!~/^non/){sum_bubble60+=1;sum_bubble_len60+=$2} else {sum_nonbubble60+=1;sum_nonbubble_len60+=$2}}
    $4>60 && $1!="total" && $4<=70 {sum70+=1;sum_len70+=$2;if($1!~/^non/){sum_bubble70+=1;sum_bubble_len70+=$2} else {sum_nonbubble70+=1;sum_nonbubble_len70+=$2}}
    $4>70 && $1!="total" && $4<=80 {sum80+=1;sum_len80+=$2;if($1!~/^non/){sum_bubble80+=1;sum_bubble_len80+=$2} else {sum_nonbubble80+=1;sum_nonbubble_len80+=$2}}
    $4>80 && $1!="total" && $4<=90 {sum90+=1;sum_len90+=$2;if($1!~/^non/){sum_bubble90+=1;sum_bubble_len90+=$2} else {sum_nonbubble90+=1;sum_nonbubble_len90+=$2}}
    $4>90 && $1!="total" && $4<=100 {sum100+=1;sum_len100+=$2;if($1!~/^non/){sum_bubble100+=1;sum_bubble_len100+=$2} else {sum_nonbubble100+=1;sum_nonbubble_len100+=$2}}
    $4>100 && $1!="total" && $4<=200 {sum200+=1;sum_len200+=$2;if($1!~/^non/){sum_bubble200+=1;sum_bubble_len200+=$2} else {sum_nonbubble200+=1;sum_nonbubble_len200+=$2}}
    $4>200{artifact+=1;artifact_len+=$2;if($1!~/^non/){artifact_bubble+=1;artifact_bubble_len+=$2} else {artifact_nonbubble+=1;artifact_nonbubble_len+=$2}} END {OFS="\t";print "0-10",sum10,sum_len10,sum_bubble10,sum_bubble_len10,sum_nonbubble10,sum_nonbubble_len10"\n""10-20",sum20,sum_len20,sum_bubble20,sum_bubble_len20,sum_nonbubble20,sum_nonbubble_len20"\n""20-30",sum30,sum_len30,sum_bubble30,sum_bubble_len30,sum_nonbubble30,sum_nonbubble_len30"\n""30-40",sum40,sum_len40,sum_bubble40,sum_bubble_len40,sum_nonbubble40,sum_nonbubble_len40"\n""40-50",sum50,sum_len50,sum_bubble50,sum_bubble_len50,sum_nonbubble50,sum_nonbubble_len50"\n""50-60",sum60,sum_len60,sum_bubble60,sum_bubble_len60,sum_nonbubble60,sum_nonbubble_len60"\n""60-70",sum70,sum_len70,sum_bubble70,sum_bubble_len70,sum_nonbubble70,sum_nonbubble_len70"\n""70-80",sum80,sum_len80,sum_bubble80,sum_bubble_len80,sum_nonbubble80,sum_nonbubble_len80"\n""80-90",sum90,sum_len90,sum_bubble90,sum_bubble_len90,sum_nonbubble90,sum_nonbubble_len90"\n""90-100",sum100,sum_len100,sum_bubble100,sum_bubble_len100,sum_nonbubble100,sum_nonbubble_len100"\n""100-200",sum200,sum_len200,sum_bubble200,sum_bubble_len200,sum_nonbubble200,sum_nonbubble_len200"\n"">200",artifact,artifact_len,artifact_bubble,artifact_bubble_len,artifact_nonbubble,artifact_nonbubble_len
    }' > result/${prefix}.${align_software}.${parameters}.md.uniq.bam.depth.mosdepth.summary.interval.txt

cat bam/${prefix}.${align_software}.${parameters}.md.uniq.bam.depth.mosdepth.summary.txt | \
    awk -F "\t" '
    $4<=50 && $1!="total" && $4>5 && $1~/^non/ {print $1}
    ' > result/${prefix}.bubblelike.list

cat bam/${prefix}.${align_software}.${parameters}.md.uniq.bam.depth.mosdepth.summary.txt | \
    awk -F "\t" '
    $4<=50 && $1!="total" && $4>5 && $1!~/^non/ {print $1}
    ' | \
    sort -t \_ -k2,2 | awk -F "\t" '{if($1~/^p/) {a[NR]=$0} else {print a[NR-1]"\t"$0}}' > result/${prefix}.bubbles.paired.list

cat bam/${prefix}.${align_software}.${parameters}.md.uniq.bam.depth.mosdepth.summary.txt | \
    awk -F "\t" '
    $4>50 && $1!="total" && $4<200 && $1~/^non/ {print $1}
    ' > result/${prefix}.real_nonbubble.list


minimap2 -cx asm20 -D -t 20 -m ${min_score} --secondary=no --cs ${INDEX} ${INDEX} > align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.paf

mawk -F "\t" 'NR==FNR {
    a[$0]=1
    next
}
    ($1 in a) && ($6 in a){
    OFS="\t"
    print $0
    }
' result/${prefix}.bubblelike.list align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.paf > align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg1.paf

mawk -F "\t" 'NR==FNR {
    a[$0]=1
    next
}
    ($1"\t"$6 in a) || ($6"\t"$1 in a){
    OFS="\t"
    print $0
    }
' result/${prefix}.bubbles.paired.list align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.paf > align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg2_bubblecontigs.paf

mawk -F "\t" 'NR==FNR {
    a[$1"\t"$6]=1
    next
}
    !($0 in a){
    OFS="\t"
    print $0
    }
' align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.paf result/${prefix}.bubbles.paired.list > result/bubble_fail.list

sort -k1,1 -k6,6 align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg2_bubblecontigs.paf > align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg2_bubblecontigs.sort.paf

python3 scripts/paf_filter_rbh.py \
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
python3 scripts/paf_filter_rbh.py \
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
    t_col=$(head -1 $tsv_file | awk -F "\t" '{for(i=1;i<=NF;i++){if($i=="t"){print i;break}}}')
    cat "$tsv_file" | sed '1d' | cut -f1,${t_col} | awk '{print $1;print $2}' | sort | uniq | awk -F '_' '{for(i=1;i<=NF;i++)if($i~/^len/){split($i,a,"n");sum+=a[2];break}}END{print sum}'
}

calculate_contig_num() {
    local tsv_file=$1
    t_col=$(head -1 $tsv_file | awk -F "\t" '{for(i=1;i<=NF;i++){if($i=="t"){print i;break}}}')
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


cat ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.tsv | sed '1d' | awk -F "\t" '{OFS="\t";print $1,$2,$3,"pair"NR;print $5,$6,$7,"pair"NR}' > ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.bed

minimap2 -cx asm20 -t 20 --secondary=no -m 200 --cs ${refrence} ${INDEX} > align/${prefix}to$(basename ${refrence} .fa).pnas.paf

cat align/${prefix}to$(basename ${refrence} .fa).pnas.paf | sort -k1,1 -k3,3n > align/${prefix}to$(basename ${refrence} .fa).pnas.sort.paf

transanno minimap2chain align/${prefix}to$(basename ${refrence} .fa).pnas.sort.paf --output align/${prefix}to$(basename ${refrence} .fa).chain

liftOver ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.bed align/${prefix}to$(basename ${refrence} .fa).chain ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.liftover.bed ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.liftover.unmap.bed

cat ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.liftover.bed | awk -F "\t" '{a[$4]+=1;if(a[$4]==2){print $4}}' > ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.liftover.list

cat ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.bed | awk -F "\t" 'NR==FNR {a[$1]=1;next;} $4 in a{OFS="\t";print;}' ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.liftover.list - | cut -f1 | sort | uniq > ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.cofirmed.list

confirmed_contig_num_bubblelikecontigs=$(cat ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.cofirmed.list | awk 'END{print NR}')
echo "Total number of confirmed bubble like contigs: $confirmed_contig_num_bubblelikecontigs"

confirmed_len_bubblelikecontigs=$(cat ${wd}/rbh_bubblelikecontigs.${min_id}.${min_cov}.cofirmed.list |  awk -F '_' '{for(i=1;i<=NF;i++)if($i~/^len/){split($i,a,"n");sum+=a[2];break}}END{print sum}')
echo "Total length of confirmed bubble like contigs: $confirmed_len_bubblelikecontigs"

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

