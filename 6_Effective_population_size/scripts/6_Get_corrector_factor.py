from Bio import SeqIO
import argparse

list_called_position=[]

with open("allchrom_all_sites_variant_unvariant_header.vcf", 'r') as f:
    for line in f:
        e=line.split("\t")
        if str(e[0][0]) != "#":
            position=e[1]
            list_called_position.append(position)

list_CDS=[]
with open("CHR*_forward_reverse_no_indiv_with_N.out", 'r') as f:
    for l in f:
        e=l.split("\t")
        contig=e[0]
        if str(contig) != "Contig_name":
            if str(contig) not in list_CDS:
                list_CDS.append(contig)

dico_CDS_callable={}
for CDS in list_CDS: # CDS_1000_chr2_244263702_244263978
    dico_CDS_callable[CDS]={}
    e=CDS.split("_")
    start=e[3] # start of a CDS
    end=e[4] # end of a CDS
    length_cds=int(end)-int(start)
    dico_CDS_callable[CDS]["length"]=length_cds
    nb_callable_sites=0
    for position in list_called_position: # if a called position is within a CDS:
        if int(position) > int(start) and int(position) < int(end): 
            nb_callable_sites+=1
        dico_CDS_callable[CDS]["nb_callable_sites"]=nb_callable_sites
        percentage_called=float(nb_callable_sites/length_cds)
        dico_CDS_callable[CDS]["precentage_called"]=percentage_called

sum_pi_N=0.0
sum_pi_S=0.0

c=0
with open("./final_all_data.txt","w") as fR:
    fR.write("Contig\tp_called\tcorr_factor\tpiN\tpiS\tpiN_corr\tpiS_corr\n")
    with open("/ALL_CHR_forward_reverse_no_indiv_with_N.out", "r") as f:
        for line in f:
            e=line.split("\t")
            contig=e[0]
            if str(e[0])!="Contig_name":
                p_called=float(dico_CDS_callable[contig]["precentage_called"])
                if float(p_called) > 0.0:
                    dico_for_bootstrap[c]={}
                    f_corr=float(1.0/float(p_called))
                    contig=e[0]
                    piN=e[26]
                    new_piN=float(piN)*float(f_corr)
                    sum_pi_N+=float(new_piN)
                    piS=e[27]
                    new_piS=float(piS)*float(f_corr)
                    sum_pi_S+=float(new_piS)
                    sum_factor+=float(p_called)
                    dico_for_bootstrap[c]["contig"]=contig
                    dico_for_bootstrap[c]["new_piN"]=new_piN
                    dico_for_bootstrap[c]["new_piS"]=new_piS
                    fR.write(str(contig)+"\t"+str(p_called)+"\t"+str(f_corr)+"\t"+str(piN)+"\t"+str(piS)+"\t"+str(new_piN)+"\t"+str(new_piS)+"\n")
                    c+=1
                        
mean_piN=float(sum_pi_N)/float(c)
mean_piS=float(sum_pi_S)/float(c)
print(mean_piN)
print(mean_piS)
ratio=float(mean_piN)/float(mean_piS)
print(ratio)

