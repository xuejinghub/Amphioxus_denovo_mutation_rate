#!/bin/bash
#SBATCH -c 20
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=inheritance
#SBATCH --output=log/inheritance_%a.out

if [ ! -e bam/Bf-${SLURM_ARRAY_TASK_ID}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam.depth.per-base.bed.gz ]
then
    cp bam/Bf-${SLURM_ARRAY_TASK_ID}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam.depth.per-base.bed.gz bam/Bf-${SLURM_ARRAY_TASK_ID}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam.depth.per-base.bed.gz
fi

python scripts/contig_postprocess_v0.9.py \
  -i result/${prefix}.contig_inheritance_v0.8.out \
  -o postprocess_result/${prefix}/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.out \
  --bubblelike test/Bf_MP_rbh_bubblelikecontigs_1to1.0.7.0.7.list \
  --len-col length \
  -r test/Bf_MP.real_nonbubble.checked.list \
  --drop_two_similar

