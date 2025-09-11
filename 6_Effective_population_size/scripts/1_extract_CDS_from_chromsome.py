from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--individual', required=True, type=str, help='Sample identifier')
args = parser.parse_args()
individual = args.individual

file_ALT = f"indi_fasta/{individual}.rename.fasta"
dico_chr_ALT={}

for record in SeqIO.parse(file_ALT, "fasta"):
    if str(record.id) not in dico_chr_ALT:
        dico_chr_ALT[str(record.id)]=record.seq

fasta_file_REVERSE=f"indi_fasta/{individual}.rename_CDS_reverse.fasta"
fasta_file_FORWARD=f"indi_fasta/{individual}.rename_CDS_forward.fasta"

GFF_CDS="bf.fmt2.gt.cds_only.SNP.allchrom.longest.gff"

# count=1
# count_reverse=0
# list_CDS=[]
# with open(fasta_file_REVERSE, "w") as fasta_to_generate_rev:
#     with open(fasta_file_FORWARD, "w") as fasta_to_generate_for:
#         with open(GFF_CDS, 'r') as gff:
#             for line in gff:
#                 e=line.split("\t")
#                 chromosome=e[0]
#                 forward_reverse=e[6]
#                 if str(chromosome) in dico_chr_ALT: 
#                     start=e[3] 
#                     end=e[4] 
#                     phase=e[7]
#                     sequence=dico_chr_ALT[str(chromosome)][int(start)-1:int(end)]  
#                     if str(forward_reverse) == "-":    
#                         if str(phase) == "0" : 
#                             sequence=sequence.reverse_complement()
#                             fasta_to_generate_rev.write(">CDS"+str(count)+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)+"\n"+str(sequence)+"\n")
#                         elif str(phase) == "1" :
#                             sequence=sequence.reverse_complement()
#                             sequence=sequence[1:]
#                             fasta_to_generate_rev.write(">CDS"+str(count)+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)+"\n"+str(sequence)+"\n")
#                         elif str(phase) == "2" : 
#                             sequence=sequence.reverse_complement()
#                             sequence=sequence[2:]
#                             fasta_to_generate_rev.write(">CDS"+str(count)+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)+"\n"+str(sequence)+"\n")
                                                          
#                     if str(forward_reverse) == "+":   
#                         if str(phase) == "0" : 
#                             fasta_to_generate_for.write(">CDS"+str(count)+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)+"\n"+str(sequence)+"\n")
#                         elif str(phase) == "1" : 
#                             sequence=sequence[1:]
#                             fasta_to_generate_for.write(">CDS"+str(count)+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)+"\n"+str(sequence)+"\n")
#                         elif str(phase) == "2" : 
#                             sequence=sequence[2:]                                
#                             fasta_to_generate_for.write(">CDS"+str(count)+"_"+str(chromosome)+"_"+str(start)+"_"+str(end)+"\n"+str(sequence)+"\n")
                                                                         
#                     count+=1
count=1
count_reverse=0
list_CDS=[]
with open(fasta_file_REVERSE, "w") as fasta_to_generate_rev:
    with open(fasta_file_FORWARD, "w") as fasta_to_generate_for:
        with open(GFF_CDS, 'r') as gff:
            # 新增字典存储基因序列片段
            gene_dict = {}
            
            for line in gff:
                e=line.split("\t")
                chromosome=e[0]
                forward_reverse=e[6]
                if str(chromosome) in dico_chr_ALT: 
                    start=int(e[3])
                    end=int(e[4])
                    phase=e[7]
                    
                    # 解析第九列获取Parent信息
                    attributes = e[8].split(';')
                    parent = [a.split('=')[1] for a in attributes if 'Parent=' in a][0]
                    
                    # 获取原始序列
                    sequence = dico_chr_ALT[str(chromosome)][start-1:end]
                    
                    # 处理相位和方向
                    if forward_reverse == "-":
                        sequence = sequence.reverse_complement()
                        if phase == "1": sequence = sequence[1:]
                        elif phase == "2": sequence = sequence[2:]
                    elif forward_reverse == "+":
                        if phase == "1": sequence = sequence[1:]
                        elif phase == "2": sequence = sequence[2:]
                    
                    # 按Parent存储序列片段
                    if parent not in gene_dict:
                        gene_dict[parent] = {
                            'strand': forward_reverse,
                            'chr': chromosome,
                            'segments': []
                        }
                    gene_dict[parent]['segments'].append((start, str(sequence)))
            
            # 处理合并后的基因序列
            for parent, data in gene_dict.items():
                # 按CDS位置排序（正链升序，负链降序）
                if data['strand'] == "+":
                    data['segments'].sort(key=lambda x: x[0])
                else:
                    data['segments'].sort(key=lambda x: -x[0])
                
                # 合并所有CDS片段
                merged_seq = ''.join([seg[1] for seg in data['segments']])
                start_pos = data['segments'][0][0]
                end_pos = data['segments'][-1][0] + len(data['segments'][-1][1]) - 1
                
                # 写入合并后的序列
                if data['strand'] == "-":
                    fasta_to_generate_rev.write(f">CDS{count}_{data['chr']}_{start_pos}_{end_pos}\n{merged_seq}\n")
                else:
                    fasta_to_generate_for.write(f">CDS{count}_{data['chr']}_{start_pos}_{end_pos}\n{merged_seq}\n")
                
                count += 1