#!/usr/bin/env python3
import pysam
import argparse
import re
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(
        description="从VCF文件中提取ALT allele的Allele Depth(AD)信息",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", required=True, help="输入VCF文件路径（需bgzip压缩和索引）")
    parser.add_argument("--output", required=True, help="输出TSV文件路径")
    args = parser.parse_args()

    try:
        vcf = pysam.VariantFile(args.input)
    except Exception as e:
        raise SystemExit(f"打开VCF文件失败: {e}")

    # 准备输出文件头
    header = ["CHROM", "POS", "QUAL", "REF", "ALT", "SAMPLE", "AD_REF", "AD_ALT", "GQ_ALT"]
    
    with open(f"{args.output}.tsv", "w") as fout:
        fout.write("\t".join(header) + "\n")
        vaf_values =[]
        for record in vcf:
            if not record.alts:
                continue
                
            # 新增染色体长度解析和位置过滤
            chrom_name = record.chrom
            match = re.search(r'_len(\d+)_', chrom_name)
            if not match:
                continue
                
            chrom_len = int(match.group(1))
            # 过滤位置在头尾200bp的变异
            if record.pos <= 200 or record.pos >= (chrom_len - 200):
                continue
            
            mq = record.info.get('MQ', '.')  # 提取MQ值
               
            # 以下保持原有样本处理逻辑不变...
            for sample_name, sample in record.samples.items():
                try:
                    ad = sample.get("AD")
                    if not ad or len(ad) < 2:
                        continue
                        
                    for alt_idx, alt in enumerate(record.alts):
                        row = [
                            record.chrom,
                            str(record.pos),
                            str(record.qual) or ".",
                            record.ref,
                            alt,
                            sample_name,
                            str(ad[0]),
                            str(ad[alt_idx+1]),
                            str(mq)
                        ]
                        if ad[alt_idx+1] < 3:
                            continue
                        fout.write("\t".join(row) + "\n")
                        vaf = ad[alt_idx+1] / (sum(ad))
                        vaf_values.append(vaf)
                        
                except KeyError:
                    continue
    # 新增可视化代码
    plt.figure(figsize=(10, 6))
    plt.hist(vaf_values, bins=50, alpha=0.7, edgecolor='black')
    plt.xlabel('Variant Allele Frequency (VAF)')
    plt.ylabel('Count')
    plt.title('VAF Distribution')
    plt.grid(True)
    plt.savefig(f"{args.output}.vaf_distribution.png", dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    main()
