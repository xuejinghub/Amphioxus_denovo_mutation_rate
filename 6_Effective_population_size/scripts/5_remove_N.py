from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--chr', required=True, type=str, help='Sample identifier')
args = parser.parse_args()
Chr = args.chr

file="results/CDS/" + str(Chr) + "_forward_reverse.fasta"
new_file="results/CDS/" + str(Chr) + "_forward_reverse_no_indiv_with_N.fasta"

with open(new_file, "w") as new:
    for record in SeqIO.parse(file, "fasta"):
        sequence=record.seq
        header=">"+str(record.id)
        if "N" not in sequence:
            new.write(header+"\n")
            new.write(str(sequence)+"\n")