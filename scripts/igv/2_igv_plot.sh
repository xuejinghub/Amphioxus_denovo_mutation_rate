#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=100G
#SBATCH -J igv_plot
#SBATCH -o log/igv_plot_%a.%j.out

prefix=Bf-${SLURM_ARRAY_TASK_ID}

singularity run -B freebayes/igv:/mnt/input -B freebayes/igv/result/:/mnt/output -B reference/:reference/ -B default/bam/:default/bam/ -B /public2/home/xuejing/bf.dnm/reassembly/read_align/bam/:/public2/home/xuejing/bf.dnm/reassembly/read_align/bam/ default/igv/igvplots_latest.sif xvfb-run --auto-servernum --server-num=1 java -jar /opt/conda/bin/igv.jar -b /mnt/input/${prefix}_igv_batch.txt