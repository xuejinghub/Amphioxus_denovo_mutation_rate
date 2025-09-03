#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=8G
#SBATCH --job-name=depth_progeny
#SBATCH --output=log/depth_progeny_%a.out

prefix=Bf-${SLURM_ARRAY_TASK_ID}

if [ ! -d result ]
then
    mkdir result
fi

if [ ! -d plot ]
then
    mkdir plot
fi

if [ ! -d tmp ]
then
    mkdir tmp
fi

zcat bam/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam.depth.per-base.bed.gz | awk -F "\t" '{if ($4 in depth) {depth[$4]+=($3-$2)} else {depth[$4]=($3-$2)}}END{for (i in depth) {print i"\t"depth[i]}}' > result/${prefix}.depth_distribution.txt

threads=$(cat plot/${prefix}.depth_distribution.out)

zcat bam/${prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam.depth.per-base.bed.gz | awk -F "\t" '{if($4 >= '$threads') {print $0}}' > tmp/${prefix}.tmp.bed

cat  ${prefix}.contig_inheritance_v0.8.postprocess.out | awk -F "\t" '{if(!($1~/^non/ && $16>0.5) && ($19~/PASS/)) {print $1}}' > tmp/${prefix}.contig_inheritance_v0.8.postprocess.list

mawk -F "\t" 'NR==FNR { 
    a[$0] = 1
    next
} 
    $1 in a { 
    print 
}' tmp/${prefix}.contig_inheritance_v0.8.postprocess.list /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/reference/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa.bed > tmp/${prefix}.contig_inheritance_v0.8.postprocess.bed

bedtools subtract -b tmp/${prefix}.tmp.bed -a tmp/${prefix}.contig_inheritance_v0.8.postprocess.bed > result/${prefix}.contig_inheritance_v0.8.postprocess.remove3foldsites.bed

cat callable_genome/progeny_depth/result/Bf-${SLURM_ARRAY_TASK_ID}.contig_inheritance_v0.8.postprocess.remove3foldsites.bed | awk -F "\t" '{sum+=($3-$2)}END{print "'${prefix}'""\t"sum}' > result/${prefix}.contig_inheritance_v0.8.postprocess.remove3foldsites.callable.stats
