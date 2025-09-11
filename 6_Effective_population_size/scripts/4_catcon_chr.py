from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--chr', required=True, type=str, help='Sample identifier')
args = parser.parse_args()
Chr = args.chr

with open("results/CDS/" + str(Chr) + "_reverse.fasta", "w") as all_cds_reverse_chr:
    for record in SeqIO.parse("results/CDS/all_CDS_reverse.fasta", "fasta"):
        header=record.id #CDS_1000_chr2_244263702_244263978|Pearl_millet|So-21-28371-02|Allele_1
        e=header.split("_")
        chromosome=e[2]
        if str(chromosome) == str(Chr):
            seq=str(record.seq)
            all_cds_reverse_chr.write(">"+header+"\n")
            all_cds_reverse_chr.write(seq+"\n") 
            
with open("results/CDS/"+ str(Chr) + "_forward.fasta", "w") as all_cds_forward_chr:
    for record in SeqIO.parse("results/CDS/all_CDS_forward.fasta", "fasta"):
        header=record.id #CDS_1000_chr2_244263702_244263978|Pearl_millet|So-21-28371-02|Allele_1
        e=header.split("_")
        chromosome=e[2]
        if str(chromosome) == str(Chr):
            seq=str(record.seq)
            all_cds_forward_chr.write(">"+header+"\n")
            all_cds_forward_chr.write(seq+"\n") 