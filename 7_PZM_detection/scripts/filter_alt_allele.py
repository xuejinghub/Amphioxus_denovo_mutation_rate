#!/usr/bin/env python3
import pysam  # type: ignore
import argparse

def check_alt_in_bam(bam, chrom, pos, alt_allele):
    """
    检查 BAM 文件中指定染色体 chrom，位置 pos（VCF 坐标，1-indexed）处是否有 read 的碱基与 alt_allele一致。
    利用 pileup 统计 0-indexed 的位点，故调用时使用 pos-1 作为起始位置。
    """
    count = 0
    for pileupcolumn in bam.pileup(chrom, pos - 1, pos, truncate=True, ignore_orphans=False, min_mapping_quality=-1 ):
        # pileupcolumn.pos 为0-indexed
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del or pileupread.is_refskip:
                    continue
                query_pos = pileupread.query_position
                if query_pos is None:
                    continue
                base = pileupread.alignment.query_sequence[query_pos]
                if base.upper() == alt_allele.upper():
                    count += 1
    return count > 1

def main():
    parser = argparse.ArgumentParser(
        description="遍历文本文件位点，检查12个BAM文件中相同位置是否存在与alt等位基因一致的read"
    )
    # 修改输入参数说明
    parser.add_argument("--in", dest="input_file", required=True,
                      help="输入文本文件路径（tab分割，列：染色体 位置 .. alt等位基因）")
    # 修改为接收12个BAM文件
    parser.add_argument("--bam", dest="bam", required=True,
                      help="12个BAM文件路径，空格分隔")
    parser.add_argument("--out", dest="out_file", required=True,  # 修改输出参数名
                      help="输出文件路径")

    args = parser.parse_args() 

    # 读取输入数据到内存
    sites = []
    with open(args.input_file, 'r') as fin:
        for idx, line in enumerate(fin):
            if idx == 0:
                header_line = line  # 保存标题行
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            sites.append({
                'chrom': parts[0],
                'pos': int(parts[1]),
                'alt': parts[4],
                'line': line  # 保留原始行内容
            })

    # 初始化结果标记（True表示需要保留）
    keep_flags = [True] * len(sites)

    # 逐个处理BAM文件节省内存
    for i in range(1,106):
        if i == 34 or i == 35 or i == int(args.bam):
            continue
        else:
            prefix = i
        bam_path = f"/public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/default/bam/Bf-{prefix}.mem2.Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.md.uniq.bam"
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for idx, site in enumerate(sites):
                if check_alt_in_bam(bam, site['chrom'], site['pos'], site['alt']):
                    keep_flags[idx] = False
    p_bam_path = "/public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/Bf_P.mem2..md.uniq.bam"
    m_bam_path = "/public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/reference/read_align/bam/Bf_M.mem2..md.uniq.bam"

    with pysam.AlignmentFile(p_bam_path, "rb") as p_bam:
        for idx, site in enumerate(sites):
            if check_alt_in_bam(p_bam, site['chrom'], site['pos'], site['alt']):
                keep_flags[idx] = False
    with pysam.AlignmentFile(m_bam_path, "rb") as m_bam:
        for idx, site in enumerate(sites):
            if check_alt_in_bam(m_bam, site['chrom'], site['pos'], site['alt']):
                keep_flags[idx] = False

    # 写入输出文件
    with open(args.out_file, 'w') as fout:
        fout.write(header_line)  # 先写入标题
        for idx, flag in enumerate(keep_flags):
            if flag:
                fout.write(sites[idx]['line'])

if __name__ == "__main__":
    main()