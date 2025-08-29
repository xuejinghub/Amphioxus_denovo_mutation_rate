#!/bin/bash
#SBATCH -c 16
#SBATCH --time=48:00:00                     
#SBATCH -J minimap2_opt
#SBATCH -o log/minimap2_opt.%J.out

threads=16
para=asm20
minimap2 \
    -ax ${para} \
    -t $threads \
    --secondary=no \
    bf.fa Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa \
    > Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.${para}.sam

samtools index Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.${para}.sort.bam

mosdepth Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.${para}.sort.bam.depth Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.${para}.sort.bam


