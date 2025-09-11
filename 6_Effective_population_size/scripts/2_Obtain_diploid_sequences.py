from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--individual', required=True, type=str, help='Sample identifier')
args = parser.parse_args()
individual = args.individual

# Forward CDS:

file_to_duplicate_seq_for=f"indi_fasta/{individual}.rename_CDS_forward.fasta"
file_duplicated_forward=f"indi_fasta/{individual}.rename_CDS_forward_DIPLOID.fasta"

no_change=['A','T','G','C']

count=1
with open(file_duplicated_forward, "w") as fR:
    for record in SeqIO.parse(str(file_to_duplicate_seq_for), "fasta"):
        sequence_1=""
        sequence_2=""
        id_seq=str(record.id) #>CDS85_chr1_3621451_3621710
        id_seq_e=id_seq.split("_")
        header=str(id_seq_e[1])+"_"+str(id_seq_e[2])+"_"+str(id_seq_e[3])
        SEQ=record.seq
        for car in SEQ:
            if str(car) == ".":
                sequence_1=str(sequence_1)+"N"
                sequence_2=str(sequence_2)+"N"
            if str(car) in no_change:
                sequence_1=str(sequence_1)+str(car)
                sequence_2=str(sequence_2)+str(car)
            if str(car) == "W":
                sequence_1=str(sequence_1)+"A"
                sequence_2=str(sequence_2)+"T"
            if str(car) == "S":
                sequence_1=str(sequence_1)+"C"
                sequence_2=str(sequence_2)+"G"
            if str(car) == "M":
                sequence_1=str(sequence_1)+"A"
                sequence_2=str(sequence_2)+"C"
            if str(car) == "K":
                sequence_1=str(sequence_1)+"G"
                sequence_2=str(sequence_2)+"T"
            if str(car) == "R":
                sequence_1=str(sequence_1)+"A"
                sequence_2=str(sequence_2)+"G"
            if str(car) == "Y":
                sequence_1=str(sequence_1)+"C"
                sequence_2=str(sequence_2)+"T"
        if  len(sequence_1) == len(sequence_2) and len(sequence_1) == len(SEQ): #check that sequences have the same length, as expected
            #>CDS_1|bf|So-21-28371-02|Allele_1
            header_1=">CDS_"+str(count)+"_"+str(header)+"|bf|"+str(individual)+"|Allele_1"
            fR.write(str(header_1)+"\n")
            fR.write(str(sequence_1)+"\n")
            header_2=">CDS_"+str(count)+"_"+str(header)+"|bf|"+str(individual)+"|Allele_2"
            fR.write(str(header_2)+"\n")
            fR.write((sequence_2)+"\n")
        else:
            print(header)
        count+=1

# Reverse CDS:

file_to_duplicate_seq_rev=f"indi_fasta/{individual}.rename_CDS_reverse.fasta"
file_duplicated_reverse=f"indi_fasta/{individual}.rename_CDS_reverse_DIPLOID.fasta"

count=1
with open(file_duplicated_reverse, "w") as fR:
    for record in SeqIO.parse(str(file_to_duplicate_seq_rev), "fasta"):
        sequence_1=""
        sequence_2=""
        id_seq=str(record.id) #>CDS85_chr1_3621451_3621710
        id_seq_e=id_seq.split("_")
        header=str(id_seq_e[1])+"_"+str(id_seq_e[2])+"_"+str(id_seq_e[3])
        SEQ=record.seq
        for car in SEQ:
            if str(car) == ".":
                sequence_1=str(sequence_1)+"N"
                sequence_2=str(sequence_2)+"N"
            if str(car) in no_change:
                sequence_1=str(sequence_1)+str(car)
                sequence_2=str(sequence_2)+str(car)
            if str(car) == "W":
                sequence_1=str(sequence_1)+"A"
                sequence_2=str(sequence_2)+"T"
            if str(car) == "S":
                sequence_1=str(sequence_1)+"C"
                sequence_2=str(sequence_2)+"G"
            if str(car) == "M":
                sequence_1=str(sequence_1)+"A"
                sequence_2=str(sequence_2)+"C"
            if str(car) == "K":
                sequence_1=str(sequence_1)+"G"
                sequence_2=str(sequence_2)+"T"
            if str(car) == "R":
                sequence_1=str(sequence_1)+"A"
                sequence_2=str(sequence_2)+"G"
            if str(car) == "Y":
                sequence_1=str(sequence_1)+"C"
                sequence_2=str(sequence_2)+"T"
        if  len(sequence_1) == len(sequence_2) and len(sequence_1) == len(SEQ): #check that sequences have the same length, as expected
            #>CDS_1|bf|So-21-28371-02|Allele_1
            header_1=">CDS_"+str(count)+"_"+str(header)+"|bf|"+str(individual)+"|Allele_1"
            fR.write(str(header_1)+"\n")
            fR.write(str(sequence_1)+"\n")
            header_2=">CDS_"+str(count)+"_"+str(header)+"|bf|"+str(individual)+"|Allele_2"
            fR.write(str(header_2)+"\n")
            fR.write((sequence_2)+"\n")
        else:
            print(header)
        count+=1