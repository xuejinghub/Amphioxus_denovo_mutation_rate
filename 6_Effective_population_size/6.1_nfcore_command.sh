
sbatch -w node80  -J batch1_bf -c 50 --mem=400G --wrap="nextflow run sarek -profile singularity --input bf_input.batch1_1.csv --fasta reference/bf.fa --igenomes_ignore --genome custom --max_memory 400.GB --tools haplotypecaller  --outdir  batch1_bf.out --max_cpus 50  --aligner bwa-mem2 --skip_tools baserecalibrator -resume disturbed_kalam --joint_germline -c sarek/conf/modules/custom.config"
