from Bio import SeqIO
import argparse
import os

# 添加命令行参数解析
parser = argparse.ArgumentParser()
parser.add_argument('--part', required=True, help='Chromosome to process')
args = parser.parse_args()

list_of_individuals=["SRR12010215_SRR12010215", "SRR12010217_SRR12010217", "SRR12010218_SRR12010218", "SRR12010219_SRR12010219", "SRR12010220_SRR12010220", "SRR12010221_SRR12010221", "SRR12010222_SRR12010222", "SRR12010223_SRR12010223", "SRR12010224_SRR12010224", "SRR12010225_SRR12010225", "SRR12010226_SRR12010226", "SRR12010228_SRR12010228", "SRR12010229_SRR12010229", "SRR12010230_SRR12010230", "SRR12010231_SRR12010231", "SRR12010232_SRR12010232", "SRR12010233_SRR12010233", "SRR12010234_SRR12010234", "SRR12010235_SRR12010235", "SRR12010236_SRR12010236", "SRR12010237_SRR12010237", "SRR12010239_SRR12010239", "SRR12010240_SRR12010240", "SRR12010241_SRR12010241", "SRR12010242_SRR12010242", "SRR12010243_SRR12010243", "SRR12010244_SRR12010244", "SRR12010246_SRR12010246", "SRR12010247_SRR12010247", "SRR12010248_SRR12010248", "SRR12010250_SRR12010250", "SRR12010251_SRR12010251", "SRR12010252_SRR12010252", "SRR12010253_SRR12010253", "SRR12010254_SRR12010254"]
# list_of_individuals=["SRR12010218_SRR12010218"]

# Getting of all the individuals for each CDS:

nb_CDS_forward=26000/2
nb_CDS_reverse=26072/2

# 创建染色体专属输出目录
output_dir = "results/CDS"
os.makedirs(output_dir, exist_ok=True)

# 修改后的处理函数
def process_cds(strand_type, start, end):
    for i in range(start, end):  # 保持原有循环范围不变
        output_file = f"{output_dir}/CDS{str(i)}_{strand_type}.fasta"
        with open(output_file, "w") as out_file:
            for indiv in list_of_individuals:
                input_file = f"indi_fasta/{indiv}.rename_CDS_{strand_type}_DIPLOID.fasta"
                for record in SeqIO.parse(input_file, "fasta"):
                    header_parts = record.id.split("_")
                    if int(header_parts[1]) == i:
                        out_file.write(f">{record.id}\n{str(record.seq)}\n")

start = (int(args.part)-1) * 500 + 1
end = int(args.part) * 500 + 1
# if int(args.part) * 500 > 13000:
#     end = 13000
# process_cds("forward", start, end)

if int(args.part) * 500 > 13036:
    end = 13036
process_cds("reverse", start, end)
