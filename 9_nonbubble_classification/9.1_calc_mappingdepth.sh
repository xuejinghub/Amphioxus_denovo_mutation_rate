#!/bin/bash
#SBATCH -w node80
#SBATCH -c 20
#SBATCH --time=48:00:00                     
#SBATCH --mem=100G
#SBATCH --job-name=align_readdepth
#SBATCH --output=log/align_readdepth_%j.out

prefix=Bf_${1}

mkdir -p wd/${prefix}

cd wd/${prefix}

if [ ! -d align ]; then
    mkdir align
fi

INDEX=${user_home}/Amphioxus/new_align/reassembly/read_align/reference/${prefix}_platanus_i3_allPhasedScaffold.rename.min150bp.fa

min_score=200

cat ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/${prefix}.mem2.${prefix}_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam.depth.mosdepth.summary.txt | \
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
    }' > ${prefix}.mem2.${prefix}_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam.depth.mosdepth.summary.interval.txt

cat ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/${prefix}.mem2.${prefix}_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam.depth.mosdepth.summary.txt | \
    awk -F "\t" '
    $4<=50 && $1!="total" && $4>5 && $1~/^non/ {print $1}
    ' > ${prefix}.bubblelike.list

cat ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/${prefix}.mem2.${prefix}_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam.depth.mosdepth.summary.txt | \
    awk -F "\t" '
    $4<=50 && $1!="total" && $4>5 && $1!~/^non/ {print $1}
    ' | \
    sort -t \_ -k2,2 | awk -F "\t" '{if($1~/^p/) {a[NR]=$0} else {print a[NR-1]"\t"$0}}' > ${prefix}.bubbles.paired.list

cat ${user_home}/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/${prefix}.mem2.${prefix}_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam.depth.mosdepth.summary.txt | \
    awk -F "\t" '
    $4>50 && $1!="total" && $4<100 && $1~/^non/ {print $1}
    ' > ${prefix}.real_nonbubble.list


minimap2 -cx asm20 -D -t 20 -m ${min_score} --secondary=no --cs ${INDEX} ${INDEX} > align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.paf

mawk -F "\t" 'NR==FNR {
    a[$0]=1
    next
}
    ($1 in a) && ($6 in a){
    OFS="\t"
    print $0
    }
' ${prefix}.bubblelike.list align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.paf > align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg1.paf

mawk -F "\t" 'NR==FNR {
    a[$0]=1
    next
}
    ($1"\t"$6 in a) || ($6"\t"$1 in a){
    OFS="\t"
    print $0
    }
' ${prefix}.bubbles.paired.list align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.paf > align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.similar_conitg2_bubblecontigs.paf

mawk -F "\t" 'NR==FNR {
    a[$1"\t"$6]=1
    next
}
    !($0 in a){
    OFS="\t"
    print $0
    }
' align/${prefix}to${prefix}.platanus_i3_allPhasedScaffold.rename.min150bp.m${min_score}.dual_yes.paf ${prefix}.bubbles.paired.list > bubble_fail.list

