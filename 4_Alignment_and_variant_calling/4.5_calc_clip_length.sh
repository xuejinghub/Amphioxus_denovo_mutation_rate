#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=100G
#SBATCH -J calc_clip_len
#SBATCH -o log/calc_clip_len_%a.%j.out

bam=bam/Bf-${SLURM_ARRAY_TASK_ID}.mem2.bf_pnas.md.uniq.bam


samtools view -F 4 ${bam} \
| gawk '{
  cigar=$6; rn=$1;
  while (match(cigar,/([0-9]+)([SH])/,a)) {
    print rn "\t" a[1] "\t" a[2];
    cigar=substr(cigar,RSTART+RLENGTH);
  }
}' > result/$(basename ${bam} .bam).clip_len.tsv

cat result/$(basename ${bam} .bam).clip_len.tsv | awk '{sum+=$2}END{print sum}' > result/$(basename ${bam} .bam).clip_len.sum

bam=bam/Bf-${SLURM_ARRAY_TASK_ID}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam

samtools view -F 4 ${bam} \
| gawk '{
  cigar=$6; rn=$1;
  while (match(cigar,/([0-9]+)([SH])/,a)) {
    print rn "\t" a[1] "\t" a[2];
    cigar=substr(cigar,RSTART+RLENGTH);
  }
}' > result/$(basename ${bam} .bam).clip_len.tsv

cat result/$(basename ${bam} .bam).clip_len.tsv | awk '{sum+=$2}END{print sum}' > result/$(basename ${bam} .bam).clip_len.sum
