#!/bin/bash
#SBATCH -w node80
#SBATCH -c 10
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=similar_align
#SBATCH --output=log/similar_align_%j.out

para=asm20
threads=10

# Process each contig in parallel array jobs
cat split/contig.list_part_${SLURM_ARRAY_TASK_ID} | while read line;
do
    parent=$(echo $line | awk -F "_" '{print $NF}')
    
    # Align modified contigs to original assembly
    minimap2 -cx ${para} -t $threads --secondary=no --cs tmp/Bf_${parent}_remove${line}_platanus_i3_allPhasedScaffold.rename.min150bp.fa  tmp/${line}.fa > tmp/${line}.${para}.paf
    
    # Sort PAF file by target contig and position
    sort -k6,6 -k8,8n tmp/${line}.${para}.paf > tmp/${line}.${para}.srt.paf
    
    # Call variants using paftools (min alignment length 100bp)
    paftools.js  call -l 100 -L 100 tmp/${line}.${para}.srt.paf > tmp/${line}.${para}.var.txt
done

# Generate alignment statistics per contig
cat split/contig.list_part_${SLURM_ARRAY_TASK_ID} | while read line;
do
    # Count unique aligned references
    grep -v "@" tmp/${line}.${para}.sam | cut -f 3 | sort | uniq | awk 'END{print "'$line'""\t" NR}' 
done > tmp/part${SLURM_ARRAY_TASK_ID}.stats

# Filter de novo mutations (DNM) against known variants
cat tmp/*.asm20.var.txt | sed '1d' | awk -F "\t" '
BEGIN {OFS="\t"}
NR==FNR {
    if ($12 == "+") {
        key = $9"\t"$11"\t"toupper($8)"\t"toupper($7)
    } else {
        c7 = $7; c8 = $8
        if (c7 == "a") c7 = "t"; else if (c7 == "t") c7 = "a";
        else if (c7 == "c") c7 = "g"; else if (c7 == "g") c7 = "c";
        if (c8 == "a") c8 = "t"; else if (c8 == "t") c8 = "a";
        else if (c8 == "c") c8 = "g"; else if (c8 == "g") c8 = "c";
        key = $9"\t"$11"\t"toupper(c8)"\t"toupper(c7)
    }
    a[key] = 1
    next
}
{
    lookup = $2"\t"$3"\t"$5"\t"$6
    if (!(lookup in a)) {
        print
    }
}' - DNM_list.removegc.tsv > test