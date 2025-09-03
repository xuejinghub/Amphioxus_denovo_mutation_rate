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
    return count > 0

def main():
    parser = argparse.ArgumentParser(
        description="遍历 VCF 位点，检查两个 BAM 文件中相同位置是否存在与 alt 等位基因一致的 read，并输出不支持 alt 的位点为 VCF 格式。"
    )
    parser.add_argument("--in", dest="vcf_file", required=True, help="输入 VCF 文件路径（必须为 bgz 压缩且建有索引）")
    parser.add_argument("--bam", dest="bam_files", action="append", required=True,
                        help="输入 BAM 文件路径，请传入两次")
    parser.add_argument("--out", dest="out_vcf", required=True, help="输出 VCF 文件路径")
    args = parser.parse_args()

    if len(args.bam_files) != 2:
        parser.error("必须提供两个 --bam 文件。")

    try:
        bam1 = pysam.AlignmentFile(args.bam_files[0], "rb")
        bam2 = pysam.AlignmentFile(args.bam_files[1], "rb")
    except Exception as e:
        parser.error(f"打开 BAM 文件失败: {e}")

    try:
        vcf_in = pysam.VariantFile(args.vcf_file)
    except Exception as e:
        parser.error(f"打开 VCF 文件失败: {e}")

    # 构建输出 VCF 文件，header 采用输入 VCF 的 header
    try:
        vcf_out = pysam.VariantFile(args.out_vcf, "w", header=vcf_in.header)
    except Exception as e:
        parser.error(f"创建输出 VCF 文件失败: {e}")

    for record in vcf_in:
        if not record.alts:
            continue

        # 这里只处理第一个 alt 等位基因；可以扩展为循环所有 alt
        alt_allele = record.alts[0]

        found_bam1 = check_alt_in_bam(bam1, record.chrom, record.pos, alt_allele)
        found_bam2 = check_alt_in_bam(bam2, record.chrom, record.pos, alt_allele)

        # 如果两个 BAM 文件均未检测到支持 alt 的 read，则将该记录写入输出 VCF
        if not found_bam1 and not found_bam2:
            vcf_out.write(record)

    bam1.close()
    bam2.close()
    vcf_in.close()
    vcf_out.close()

if __name__ == "__main__":
    main()