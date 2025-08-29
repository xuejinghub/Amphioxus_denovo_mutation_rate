#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=100G
#SBATCH -J inhe_plot
#SBATCH -o log/inhe_lpt_%a.%j.out

prefix=Bf-${SLURM_ARRAY_TASK_ID}

if [ ! -d /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/cosin_similarity/postprocess ]
then
    mkdir /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/cosin_similarity/postprocess
fi

python /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/cosin_similarity/contig_postprocess.py \
    --in /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/cosin_similarity/result/${prefix}.contig_inheritance_v0.8.out \
    --out /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/cosin_similarity/postprocess/${prefix}.contig_inheritance_v0.8.postprocess.out

cat /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/cosin_similarity/postprocess/${prefix}.contig_inheritance_v0.8.postprocess.out | grep "PASS" | cut -f1  > pass_list/${prefix}.passlist

mawk -F "\t" 'NR==FNR { 
    a[$0] = 1
    next
} 
    $5 in a { 
    print 
}' pass_list/${prefix}.passlist Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.site_res.liftover.merged2.bed > tmp/${prefix}.tmp.list.pass

Rscript /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/default/chromosome_plot/plot_parental_source.r -i /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/default/inheritance_filter/tmp/${prefix}.tmp.list.pass -o /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/default/chromosome_plot/plot/${prefix}.pdf -p ${prefix}