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

python contig_analysis_v0.8.py -p Bf_P.mem2..md.uniq.bam.depth.per-base.bed.gz  -m Bf_M.mem2..md.uniq.bam.depth.per-base.bed.gz  --progeny bam/Bf-${SLURM_ARRAY_TASK_ID}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam.depth.per-base.bed.gz  --threshold 0.9 --pr_id Bf-${SLURM_ARRAY_TASK_ID} --processes 20 -o result/Bf-${SLURM_ARRAY_TASK_ID}.contig_inheritance_v0.8.out --method hybrid

