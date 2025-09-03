#!/bin/bash
#SBATCH -c 5
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=similar_align
#SBATCH --output=log/similar_align_%j.out

if [ ! -d align ]
then
    mkdir align
fi

parent=$1
if [ "$parent" == "M" ]; then
    opposite="P"
else
    opposite="M"
fi

# Cross-align parental genomes using minimap2
minimap2 -cx asm20 -t 50 --secondary=no --cs Bf_M_platanus_i3_allPhasedScaffold.rename.min150bp.fa Bf_P_platanus_i3_allPhasedScaffold.rename.min150bp.fa  > align/Bf-PtoM.platanus_i3_allPhasedScaffold.rename.min150bp.paf
minimap2 -cx asm20 -t 50 --secondary=no --cs Bf_P_platanus_i3_allPhasedScaffold.rename.min150bp.fa Bf_M_platanus_i3_allPhasedScaffold.rename.min150bp.fa  > align/Bf-MtoP.platanus_i3_allPhasedScaffold.rename.min150bp.paf

# Filter PAF file using pre-defined similar contig list
mawk -F "\t" 'NR==FNR {
    a[$0]=1
    next
}
    ($1 in a) && ($6 in a){
    OFS="\t"
    print $0
    }
' similar_conitg.list align/Bf-${parent}to${opposite}.platanus_i3_allPhasedScaffold.rename.min150bp.paf > align/Bf-${parent}to${opposite}.platanus_i3_allPhasedScaffold.rename.min150bp.similar_conitg1.paf

# Process contig pairs in parallel array jobs
cat split/similar_conitg.list_part${SLURM_ARRAY_TASK_ID} | while read contig; do 
    parent=$(echo $contig | mawk '{if($0~/M$/) {print "M"} else if($0~/P$/){print "P"}}')
    if [ $parent == "M" ]
    then 
        parent_other="P"
    else 
        parent_other="M" 
    fi
    # Find reciprocal best matches between parental contigs
    mawk -F "\t" -v contig="${contig}" '{if(contig==$1)print $6}' align/Bf-${parent}to${parent_other}.platanus_i3_allPhasedScaffold.rename.min150bp.similar_conitg.paf | while read opposite; do
        mawk -F "\t" -v opposite="${opposite}" '{if(opposite==$1)print $6}' align/Bf-${parent_other}to${parent}.platanus_i3_allPhasedScaffold.rename.min150bp.similar_conitg.paf | while read other; do
            mawk -F "\t" -v contig="${contig}" -v other="${other}" -v opposite="${opposite}" -v parent="${parent}" 'BEGIN{if(contig==other){if(parent=="M") {print contig"\t"opposite} else if (parent=="P") {print opposite"\t"contig}}}'
        done
    done
done > tmp/similar_conitg_part${SLURM_ARRAY_TASK_ID}.out

# Generate mapping quality statistics using samtools
samtools view  -O SAM bam/Bf_${parent}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam | grep -v "@" | mawk -F "\t" '{if($3 in sum) {sum[$3]+=$5;count[$3]+=1} else {sum[$3]=$5;count[$3]=1}}END{for(i in sum){print i"\t"sum[i]"\t"count[i]"\t"sum[i]/count[i]}}' > Bf_${parent}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam.mapq
samtools view  -O SAM bam/Bf_${parent}.mem2.Bf_${parent}_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam | grep -v "@" | mawk -F "\t" '{if($3 in sum) {sum[$3]+=$5;count[$3]+=1} else {sum[$3]=$5;count[$3]=1}}END{for(i in sum){print i"\t"sum[i]"\t"count[i]"\t"sum[i]/count[i]}}' > Bf_${parent}.mem2.Bf_${parent}_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam.mapq