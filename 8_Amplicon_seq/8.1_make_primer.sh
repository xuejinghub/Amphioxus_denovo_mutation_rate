#!/bin/bash
#SBATCH -c 10
#SBATCH --time=48:00:00
#SBATCH --mem=20G
#SBATCH --job-name=make_primer
#SBATCH --output=log/make_primer_%a.out

length=$1

wd=validate/run250614

if [ ! -d callable_genome/individual ]
then
    mkdir -p callable_genome/individual
fi

prefix=Bf-${SLURM_ARRAY_TASK_ID}

if [ ! -d callable_genome/individual/${prefix} ]
then
    mkdir -p callable_genome/individual/${prefix}  
fi

if [ ! -d primer/${prefix}/${length} ]
then
    mkdir -p primer/${prefix}/${length}
fi

if [ ! -d pcr_products/${prefix}/${length}  ]
then
    mkdir -p pcr_products/${prefix}/${length}
fi

## Print basic information in log:

echo -e "Working dictory: $wd\nSample ID: ${prefix}\nAmplicon Length: ${length}\nCPUs used: $SLURM_CPUS_ON_NODE\n"

if [ ! -e callable_genome/individual/${prefix}/${prefix}.del_homo_contig.list ]
then
    cat ${wd}/bf_germlineDNM_seleceted_somaticDNM.list.tsv | awk -F "\t" '$1=="'${SLURM_ARRAY_TASK_ID}'"' | while read line
    do
    contig=$(echo $line | sed 's/ /\t/g' | awk -F "\t" '{print $2}')
    cat validate/run250613/lo_ref_result.germlineAndsomaticDNM.tsv | grep "${contig}" | awk -F "\t" '{for(i=4;i<=7;i++){if(($i!~/^NA/) && ($i!~/^'$contig'/)) {print $i}}}'
    done | cut -d \: -f1 | grep -v -f add_primers.contiglist.tsv > callable_genome/individual/${prefix}/${prefix}.del_homo_contig.list
fi

# Optimize for clusteredDNM
if [ ! -e primer/${prefix}/${prefix}.primer_regions.txt ]
then
    cat ${wd}/bf_germlineDNM_seleceted_somaticDNM.list.tsv ${wd}/add_primers.sitelist.fmt.tsv | awk -F "\t" '$1=="'${SLURM_ARRAY_TASK_ID}'"' | awk -F "\t" '{OFS="\t";if($4=="Cluster"){if($2 in contig) {contig[$2]+=$3;contig_num[$2]+=1} else {contig[$2]=$3;contig_num[$2]=1}} else if($4=="PASS") {print $2,$3,$3+1,$2":"$3,"","B","NotW"}}END{for(i in contig) {pos=contig[i]/contig_num[i]; print i,int(pos),int(pos)+1,i":"int(pos),"","B","NotW"}}' | sort -k1,1 > primer/${prefix}/${prefix}.primer_regions.txt
fi

if [ ! -e callable_genome/individual/${prefix}/${prefix}.contig_inherited.fa ]
then
    seqkit grep -v -f callable_genome/individual/${prefix}/${prefix}.del_homo_contig.list reference/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa > callable_genome/individual/${prefix}/${prefix}.contig_inherited.fa
fi

if [ ! -e callable_genome/individual/${prefix}/${prefix}.contig_inherited.fa.fai ]
then
    bwa index callable_genome/individual/${prefix}/${prefix}.contig_inherited.fa
    samtools faidx callable_genome/individual/${prefix}/${prefix}.contig_inherited.fa
fi


cd primer/${prefix}/${length}

cp ../${prefix}.primer_regions.txt .

python validate/NGS-PrimerPlex/NGS_primerplex.py --regions-file ${prefix}.primer_regions.txt --reference-genome ${wd}/callable_genome/individual/${prefix}/${prefix}.contig_inherited.fa --min-amplicon-length 100 --max-amplicon-length 300 --optimal-amplicon-length ${length} --min-primer-length 18 --max-primer-length 30 --min-primer-melting-temp 60 --max-primer-melting-temp 68 --optimal-primer-melting-temp 64 --min-primer-gc 30 --max-primer-gc 70 --optimal-primer-gc 50 --primers-number1 50 --threads 10 --return-variants-number 200

for i in {1..200}
do
    python xlsx_to_txt.py ${prefix}.primer_regions_primers_combination_${i}_info.xls ${prefix}.primer_regions_primers_combination_${i}_info.txt
done

cat ${wd}/bf_germlineDNM_seleceted_somaticDNM.list.tsv ${wd}/add_primers.sitelist.fmt.tsv | awk -F "\t" '$1=="'${SLURM_ARRAY_TASK_ID}'"' | while read line
do
    contig=$(echo $line | sed 's/ /\t/g' | awk -F "\t" '{print $2}')
    cat validate/run250613/lo_ref_result.germlineAndsomaticDNM.tsv | grep "${contig}" | awk -F "\t" '{for(i=4;i<=7;i++){if($i!~/^NA/) {print $i}}}'| cut -d \: -f1 > ${wd}/callable_genome/individual/${prefix}/${prefix}.${contig}_homolog.list
    if [ -s ${wd}/callable_genome/individual/${prefix}/${prefix}.${contig}_homolog.list ]
    then 
        if [ ! -e ${wd}/callable_genome/individual/${prefix}/${prefix}.${contig}_homolog.fa ]
        then
            seqkit grep -f ${wd}/callable_genome/individual/${prefix}/${prefix}.${contig}_homolog.list reference/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa > ${wd}/callable_genome/individual/${prefix}/${prefix}.${contig}_homolog.fa
        fi
        for i in {1..200}
        do
            primer1=$(cat ${prefix}.primer_regions_primers_combination_${i}_info.txt | grep "${contig}" | awk -F "\t" '{print $2}')
            primer2=$(cat ${prefix}.primer_regions_primers_combination_${i}_info.txt | grep "${contig}" | awk -F "\t" '{print $3}')
            [ -z "${primer1}" ] && continue
            echo -e "${primer1}\t${primer2}\t${prefix}_${contig}"        
        done | sort | uniq | awk '{OFS="\t";print $0"_"NR}'> ${wd}/pcr_products/${prefix}/${length}/${prefix}.${contig}.len${length}.pairs 
        perl ${wd}/in_silico_PCR.pl -s ${wd}/callable_genome/individual/${prefix}/${prefix}.${contig}_homolog.fa -p ${wd}/pcr_products/${prefix}/${length}/${prefix}.${contig}.len${length}.pairs -m 2 -i -f ${wd}/pcr_products/${prefix}/${length}/${prefix}.${contig}.len${length}.product.fa > ${wd}/pcr_products/${prefix}/${length}/${prefix}.${contig}.len${length}.out
    fi
done 

cat ${wd}/bf_germlineDNM_seleceted_somaticDNM.list.tsv ${wd}/add_primers.sitelist.fmt.tsv | awk -F "\t" '$1=="'${SLURM_ARRAY_TASK_ID}'"' | while read line
do
    contig=$(echo $line | sed 's/ /\t/g' | awk -F "\t" '{print $2}')
    for i in 1
    do
        primer1=$(cat ${prefix}.primer_regions_primers_combination_${i}_info.txt | grep "${contig}" | awk -F "\t" '{print $2}')
        primer2=$(cat ${prefix}.primer_regions_primers_combination_${i}_info.txt | grep "${contig}" | awk -F "\t" '{print $3}')
        [ -z "${primer1}" ] && continue
        echo -e "${primer1}\t${primer2}\t${prefix}_${contig}"        
    done | sort | uniq | awk '{OFS="\t";print $0"_"NR}'> ${wd}/pcr_products/${prefix}/${length}/${prefix}.${contig}.len${length}.pair1.pairs 
    perl ${wd}/in_silico_PCR.pl -s ${wd}/callable_genome/individual/${prefix}/${prefix}.${contig}_homolog.fa -p ${wd}/pcr_products/${prefix}/${length}/${prefix}.${contig}.len${length}.pair1.pairs -m 2 -i -f ${wd}/pcr_products/${prefix}/${length}/${prefix}.${contig}.len${length}.product.pair1.fa > ${wd}/pcr_products/${prefix}/${length}/${prefix}.${contig}.len${length}.pair1.out
done 


cat ${wd}/bf_germlineDNM_seleceted_somaticDNM.list.tsv | awk -F "\t" '$1=="'${SLURM_ARRAY_TASK_ID}'"' | while read line
do
    contig=$(echo $line | sed 's/ /\t/g' | awk -F "\t" '{print $2}')
    perl ${wd}/in_silico_PCR.pl -s reference/Bf_MP_platanus_i3_allPhasedScaffold.rename.min150bp.fa -p ${wd}/pcr_products/${prefix}/${length}/${prefix}.${contig}.len${length}.pairs  -m 2 -i -f ${wd}/pcr_products/${prefix}/${prefix}.${contig}.all_genome.product.fa > ${wd}/pcr_products/${prefix}/${prefix}.${contig}.all_genome.out
done 