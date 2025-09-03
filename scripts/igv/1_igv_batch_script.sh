#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=100G
#SBATCH -J igv_plot
#SBATCH -o log/igv_batch_script_%a.%j.out
# 1. 定义变量（根据实际情况修改）
prefix=Bf-${SLURM_ARRAY_TASK_ID}
if [ ! -d result/${prefix} ]
then
    mkdir -p result/${prefix}
fi
SNAP_DIR=result/${prefix}

# 2. 创建批处理文件头部（加载文件）
echo -e "new" > Bf-${SLURM_ARRAY_TASK_ID}_igv_batch.txt
echo -e "genome reference/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa" >> Bf-${SLURM_ARRAY_TASK_ID}_igv_batch.txt
echo -e "load default/bam/Bf_P.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam" >> Bf-${SLURM_ARRAY_TASK_ID}_igv_batch.txt
echo -e "load default/bam/Bf_M.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam" >> Bf-${SLURM_ARRAY_TASK_ID}_igv_batch.txt
for i in {0..9}
do
    id=Bf-$[SLURM_ARRAY_TASK_ID + i]
    echo -e "load default/bam/${id}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam" >> Bf-${SLURM_ARRAY_TASK_ID}_igv_batch.txt
done
echo -e "snapshotDirectory $SNAP_DIR" >> Bf-${SLURM_ARRAY_TASK_ID}_igv_batch.txt

# 3. 处理VCF文件并追加命令到批处理文件
for i in {0..9}
do
    id=Bf-$[SLURM_ARRAY_TASK_ID + i]
    # VCF_FILE=postzygotic/filter/${id}.alt_ad.filterrepeat.filteroverlap.tsv            # VCF文件名
    VCF_FILE=gc/filter_gc/result/${id}.mem2.Bf_MP_platanus_i3_allPhasedScaffold_phased.md.uniq.freebayes.noMP.filtered.removegc2.vcf          # VCF文件名
    cat $VCF_FILE | head -1 | awk '!/^#/ {
        chr = $1; 
        pos = $2; 
        start = (pos - 20 < 1) ? 1 : pos - 20; 
        end = pos + 20; 
        print "goto " chr ":" start "-" end; 
        print "snapshot " "'${id}'_" chr "_" pos "_region.png"; 
        print "goto " chr; 
        print "sleep 1"  # 新增等待确保窗口调整完成 
        print "snapshot " "'${id}'_" chr "_" pos "_zoomed.png"
    }'  >> Bf-${SLURM_ARRAY_TASK_ID}_igv_batch.txt
done
# 4. 追加批处理文件的尾部（退出）
echo -e "exit" >> Bf-${SLURM_ARRAY_TASK_ID}_igv_batch.txt
