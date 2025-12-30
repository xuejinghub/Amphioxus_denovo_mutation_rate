#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=100G
#SBATCH --job-name=1generate_sim_reads
#SBATCH --output=log/1generate_sim_reads_%j.out

prefix=Bf-${SLURM_ARRAY_TASK_ID}
refrence=${user_home}/Amphioxus/bf.dnm/reference/bf.fa
parents_contigs=${user_home}/Amphioxus/new_align/reassembly/read_align/reference/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa
mut=1000
depth=40
CG_dir=${user_home}/Amphioxus/new_align/reassembly/read_align/non_bubble/inheritance/similarity_analysis/postprocess_result/${prefix}
fastq_dir=simulation_data
uncallableE3_pair=${CG_dir}/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.uncallableE3_pairs.tsv
uncallableE3_callable=${CG_dir}/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.uncallableE3.remove3foldsites.bed

if [ ! -d sim_fasta/${prefix} ]; then
    mkdir -p sim_fasta/${prefix}
fi
if [ ! -e Bf_M.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.site.bed ]; then
zcat ${user_home}/Amphioxus/new_align/reassembly/read_align/freebayes_update1/vcf/Bf_M.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz | grep -v "#" | mawk -F "\t" '{OFS="\t";print $1,$2-1,$2}' | sort -k1,1 -k2,2n > ${user_home}/Amphioxus/new_align/reassembly/read_align/FNR/Bf_M.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.site.bed
fi
if [ ! -e Bf_P.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.site.bed ]; then
zcat ${user_home}/Amphioxus/new_align/reassembly/read_align/freebayes_update1/vcf/Bf_P.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.vcf.gz | grep -v "#" | mawk -F "\t" '{OFS="\t";print $1,$2-1,$2}' | sort -k1,1 -k2,2n > ${user_home}/Amphioxus/new_align/reassembly/read_align/FNR/Bf_P.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.site.bed
fi

if [ ! -e sim_fasta/${prefix}/${prefix}.inherited_contigs.single_contigE3.sim_mut${mut}.fa ]; then
    cat ${uncallableE3_pair} | awk -F "\t" '{split($1,M,"_");for(i=1;i<=length(M);i++)if(M[i]~/^len/){split(M[i],a,"n")};split($2,P,"_");for(i=1;i<=length(P);i++)if(P[i]~/^len/){split(P[i],b,"n")};if(a[2]>b[2]){print $1} else if(a[2]<b[2]){print $2}}' | sort | uniq > callable_genome/bed_update/$(basename ${uncallableE3_pair} _pairs.tsv).single_contig.list
    cat ${uncallableE3_callable}  | awk -F "\t" 'NR==FNR{a[$1]=1;next;} $1 in a{OFS="\t";print;}' callable_genome/bed_update/$(basename ${uncallableE3_pair} _pairs.tsv).single_contig.list - > callable_genome/bed_update/$(basename ${uncallableE3_callable} .tsv).single_contig.bed
    cut -f1 ${CG_dir}/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.callable_filtered.remove3foldsites.remove_unknow_phased_contig.bed | sort | uniq  > callable_genome/bed_update/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.callable_filtered.contig.list
    cat callable_genome/bed_update/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.callable_filtered.contig.list callable_genome/bed_update/$(basename ${uncallableE3_pair} _pairs.tsv).single_contig.list | sort | uniq > callable_genome/bed_update/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.all.contig.list
    seqkit grep -f callable_genome/bed_update/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.all.contig.list ${parents_contigs} > sim_fasta/${prefix}/${prefix}.inherited_contigs.single_contigE3.fa
    cat ${CG_dir}/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.callable_filtered.remove3foldsites.remove_unknow_phased_contig.bed callable_genome/bed_update/$(basename ${uncallableE3_callable} .tsv).single_contig.bed | sort -k1,1 -k2,2n | bedtools subtract -a - -b Bf_P.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.site.bed Bf_M.mem2.Bf_MP_platanus_i3_allPhasedScaffold.md.uniq.freebayes.site.bed -sorted > callable_genome/bed_update/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.callable_filtered.uncallableE3.remove3foldsites.removeParentsVariants.bed
    python replace_v2.0.py \
        sim_fasta/${prefix}/${prefix}.inherited_contigs.single_contigE3.fa \
        callable_genome/bed_update/${prefix}.contig_inheritance_v0.8.postprocess_v0.9.callablegenome_v0.3.callable_filtered.uncallableE3.remove3foldsites.removeParentsVariants.bed \
        ${mut} \
        -p sim_fasta/${prefix}/${prefix}.inherited_contigs.single_contigE3.sim_mut${mut}
fi

source /public/home/xuejing/mambaforge/etc/profile.d/conda.sh
conda activate art

if [ ! -e $fastq_dir/${prefix}_1.fq ] || [ ! -e $fastq_dir/${prefix}_2.fq ]
then
art_illumina -p -m 250 -s 50 -f $depth -l 150 -ef -ss HS25 -i sim_fasta/${prefix}/${prefix}.inherited_contigs.single_contigE3.sim_mut${mut}.fa -o $fastq_dir/${prefix}_mut${mut}_
fi

if [ ! -e $fastq_dir/${prefix}_1.fq.gz ] || [ ! -e $fastq_dir/${prefix}_2.fq.gz ]
then
    gzip $fastq_dir/${prefix}_mut${mut}_1.fq
    gzip $fastq_dir/${prefix}_mut${mut}_2.fq
fi





