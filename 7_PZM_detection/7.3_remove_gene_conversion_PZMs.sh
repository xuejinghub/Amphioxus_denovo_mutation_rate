#!/bin/bash
#SBATCH -w node80
#SBATCH -c 10
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=self_align
#SBATCH --output=log/self_align_%j.out

para=asm20
threads=10
wd=pzDNM_update

if [ ! -d $wd ]
then
    mkdir $wd   
fi

cat ${wd}/split/contig.list_pzDNM_part_${SLURM_ARRAY_TASK_ID} | cut -f 1 | while read line;
do
    parent=$(echo $line | awk -F "_" '{print $NF}')
    seqkit grep -p $line /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/reference/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa > ${wd}/${line}.fa
    seqkit grep -v -p $line Bf_${parent}_platanus_i3_allPhasedScaffold.rename.min150bp.fa > $wd/Bf_${parent}_remove${line}_platanus_i3_allPhasedScaffold.rename.min150bp.fa
    minimap2 -cx ${para} -t $threads --secondary=no --cs ${wd}/Bf_${parent}_remove${line}_platanus_i3_allPhasedScaffold.rename.min150bp.fa  ${wd}/${line}.fa > ${wd}/${line}.${para}.paf
    rm ${wd}/Bf_${parent}_remove${line}_platanus_i3_allPhasedScaffold.rename.min150bp.fa
    sort -k6,6 -k8,8n ${wd}/${line}.${para}.paf > ${wd}/${line}.${para}.srt.paf
    paftools.js  call -l 100 -L 100 ${wd}/${line}.${para}.srt.paf > ${wd}/${line}.${para}.var.txt
done

cat split/contig.list_pzDNM_part_${SLURM_ARRAY_TASK_ID} | while read line;
do
    grep -v "@" ${wd}/${line}.${para}.sam | cut -f 3 | sort | uniq | awk 'END{print "'$line'""\t" NR}' 
done > ${wd}/part${SLURM_ARRAY_TASK_ID}.stats

cat ${wd}/*.asm20.var.txt | awk -F "\t" '$1=="V"'|  awk -F "\t" '
BEGIN {OFS="\t"}
NR==FNR {
    if ($12 == "+") {
        key = $9"\t"$11"\t"toupper($8)"\t"toupper($7)
    } else {
        c7 = $7; c8 = $8
        if (c7 == "a") c7 = "t"; else if (c7 == "t") c7 = "a";
        else if (c7 == "c") c7 = "g"; else if (c7 == "g") c7 = "c";
        if (c8 == "a") c8 = "t"; else if (c8 == "t") c8 = "a";
        else if (c8 == "c") c8 = "g"; else if (c8 == "g") c8 = "c";
        key = $9"\t"$11"\t"toupper(c8)"\t"toupper(c7)
    }
    a[key] = 1
    next
}
FNR == 1 { next }
{
    lookup = $1"\t"$2"\t"$4"\t"$5
    if (!(lookup in a)) {
        print
    }
}' -  postzygoticDNM_removeSimilarnonbubblehighDepth.info.tsv > ${wd}/postzygoticDNM_removeSimilarnonbubblehighDepth.info.result.removegc.tsv

cat ${wd}/*.asm20.var.txt pzDNM/*.asm20.var.txt | awk -F "\t" '$1=="V"'|  awk -F "\t" '
BEGIN {OFS="\t"}
NR==FNR {
    if ($12 == "+") {
        key = $9"\t"$11"\t"toupper($8)"\t"toupper($7)
    } else {
        c7 = $7; c8 = $8
        if (c7 == "a") c7 = "t"; else if (c7 == "t") c7 = "a";
        else if (c7 == "c") c7 = "g"; else if (c7 == "g") c7 = "c";
        if (c8 == "a") c8 = "t"; else if (c8 == "t") c8 = "a";
        else if (c8 == "c") c8 = "g"; else if (c8 == "g") c8 = "c";
        key = $9"\t"$11"\t"toupper(c8)"\t"toupper(c7)
    }
    a[key] = 1
    next
}
FNR == 1 { next }
{
    lookup = $1"\t"$2"\t"$4"\t"$5
    if (!(lookup in a)) {
        print
    }
}' -  All_samples.alt_ad.filterrepeat.filteroverlap.filteralt_allsamples.tsv > All_samples.alt_ad.filterrepeat.filteroverlap.filteralt_allsamples.removegc.tsv