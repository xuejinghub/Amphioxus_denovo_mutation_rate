#!/bin/bash
#SBATCH -w node79
#SBATCH -c 10
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=2inheritance_simulation
#SBATCH --output=log/2inheritance_simulation_%j_%a.out

chr=Chr3
het=0.03
depth=40
species=bf
batch=${1}
prefix=${chr}_sim_het${het}_F1batch${batch}_dep${depth}
genome_size=$(cat ${user_home}/Amphioxus/bf.dnm/reference/bf.fa.fai | awk -v chr="${chr}" '$1==chr{print $2}')
refrence=${user_home}/Amphioxus/bf.dnm/reference/bf.fa
cpu=20
fastq_dictory=${user_home}/Amphioxus/new_align/reassembly/read_align/cosin_similarity/simulation/wd/${prefix}/fq
min_score=200
align_software=mem2
parameters=${chr}_sim_het${het}_MP_dep${depth}_allPhasedScaffold.min150.rename
threads=${cpu}
script_dir=${user_home}/Amphioxus/new_align/reassembly/read_align/cosin_similarity/simulation/script

sim_genome_dir=${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/align_readdepth/sim/pisson_assebly/sim_genome

if [ ! -d wd/${prefix} ]; then
    mkdir -p wd/${prefix}
fi

cd wd/${prefix}

if [ ! -d tmp ]; then
    mkdir tmp
fi

if [ ! -d fq ]; then
    mkdir fq
fi

if [ ! -d fa ]; then
    mkdir fa
fi

if [ ! -d bam ]; then
    mkdir bam
fi

if [ ! -d align ]; then
    mkdir align
fi

if [ ! -d result ]; then
    mkdir result
fi

cat ${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/align_readdepth/sim/pisson_assebly/wd/Chr3_sim_het0.03_P_dep40_batch1/Chr3_sim_het0.03_P_dep40_batch1.real_nonbubble.checked.list | awk '{print $1"_Bf_P"}' > Chr3_sim_het0.03_P_dep40_batch1.real_nonbubble.checked.list

cat ${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/align_readdepth/sim/pisson_assebly/wd/Chr3_sim_het0.03_M_dep40_batch2/Chr3_sim_het0.03_M_dep40_batch2.real_nonbubble.checked.list | awk '{print $1"_Bf_M"}' > Chr3_sim_het0.03_M_dep40_batch2.real_nonbubble.checked.list

cat Chr3_sim_het0.03_P_dep40_batch1.real_nonbubble.checked.list Chr3_sim_het0.03_M_dep40_batch2.real_nonbubble.checked.list > Chr3_sim_het0.03_MP_dep${depth}.real_nonbubble.checked.list

cut -f1 ${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/align_readdepth/sim/pisson_assebly/wd/Chr3_sim_het0.03_M_dep40_batch2/result/min_id0.7_min_cov0.7/rbh_bubblelikecontigs.0.7.0.7.bed | sort | uniq | awk '{print $1"_Bf_M"}' > Chr3_sim_het0.03_M_dep40_batch2_rbh_bubblelikecontigs.0.7.0.7.list

cut -f1 ${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/align_readdepth/sim/pisson_assebly/wd/Chr3_sim_het0.03_P_dep40_batch1/result/min_id0.7_min_cov0.7/rbh_bubblelikecontigs.0.7.0.7.bed | sort | uniq | awk '{print $1"_Bf_M"}' > Chr3_sim_het0.03_P_dep40_batch2_rbh_bubblelikecontigs.0.7.0.7.list

cat Chr3_sim_het0.03_M_dep40_batch2_rbh_bubblelikecontigs.0.7.0.7.list Chr3_sim_het0.03_P_dep40_batch2_rbh_bubblelikecontigs.0.7.0.7.list > Chr3_sim_het0.03_MP_dep${depth}_rbh_bubblelikecontigs.0.7.0.7.list

minimap2 -cx asm20 -t 50 --secondary=no --cs fa/${chr}_sim_het${het}_P_dep40_allPhasedScaffold.min150.rename.fa fa/${chr}_sim_het${het}_M_dep40_allPhasedScaffold.min150.rename.fa > align/${chr}_sim_het${het}_MtoP_dep${depth}_allPhasedScaffold.min150.rename.paf

minimap2 -cx asm20 -t 50 --secondary=no --cs fa/${chr}_sim_het${het}_M_dep40_allPhasedScaffold.min150.rename.fa fa/${chr}_sim_het${het}_P_dep40_allPhasedScaffold.min150.rename.fa > align/${chr}_sim_het${het}_PtoM_dep${depth}_allPhasedScaffold.min150.rename.paf

python $script_dir/rbh_blocks.py align/${chr}_sim_het${het}_MtoP_dep${depth}_allPhasedScaffold.min150.rename.paf align/${chr}_sim_het${het}_PtoM_dep${depth}_allPhasedScaffold.min150.rename.paf  --out-prefix mat_pat_rbh --min-ident-block 0.7 --overlap-frac 0.5

cat mat_pat_rbh.contig_pairs.tsv | cut -f1,3 > mat_pat_rbh.contig_pairs.list

python $script_dir/contig_analysis_v0.8.py -p bam/Chr3_sim_het0.03_P_dep40_batch1.${align_software}.${parameters}.md.uniq.bam.depth.per-base.bed.gz  -m bam/Chr3_sim_het0.03_M_dep40_batch2.${align_software}.${parameters}.md.uniq.bam.depth.per-base.bed.gz  --progeny bam/${prefix}.${align_software}.${parameters}.md.uniq.bam.depth.per-base.bed.gz  --threshold 0.5 --pr_id ${batch} --processes 10 -o result/${prefix}.contig_inheritance_v0.8.out --method hybrid


if [ ! -d postprocess_result/ ]; then
  mkdir -p postprocess_result/
fi

python $script_dir/contig_postprocess_v0.9.py \
  -i result/${prefix}.contig_inheritance_v0.8.out \
  -o postprocess_result/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.out \
  --bubblelike Chr3_sim_het0.03_MP_dep${depth}_rbh_bubblelikecontigs.0.7.0.7.list \
  --len-col length \
  -r Chr3_sim_het0.03_MP_dep${depth}.real_nonbubble.checked.list \
  --drop_two_similar
  
python $script_dir/callable_genome.v0.3.py \
  -i postprocess_result/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.out \
  -o postprocess_result/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3 \
  --parental-homology mat_pat_rbh.contig_pairs.list --len-col length 

bedtools subtract -a postprocess_result/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.callable_filtered.bed \
    -b ${user_home}/Amphioxus/new_align/reassembly/read_align/default/callable_genome/progeny_depth/tmp/${prefix}.tmp.bed > postprocess_result/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.callable_filtered.remove3foldsites.bed

