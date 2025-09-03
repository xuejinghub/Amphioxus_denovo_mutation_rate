#!/bin/bash
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=similar_count
#SBATCH --output=log/similar_count_%a.out


prefix=Bf-${SLURM_ARRAY_TASK_ID}
cat ${prefix}.contig_inheritance_v0.8.postprocess.out | awk -F '\t' '
{
    match($1, /bubble[0-9]+/);
    bubble_num = substr($1, RSTART+6, RLENGTH-6);
    parent = substr($1, length($1));
    
    key = parent "_" bubble_num;
    
    if ($1 ~ /^primary/) {
        primary[key] = $16;
        primary_line[key] = $0;
        primary_depth[key] = $19;
        primary_len[key] = $17
    }
    else if ($1 ~ /^secondary/) {
        if (key in primary && 
            primary[key] > 0.25 && $16 > 0.25 &&
            primary_depth[key] != "low_depth" && $19 != "low_depth") {
            print primary_line[key];
            print $0;
            count++;
            if($19=="PASS"){len+=$17;} else {len+=primary_len[key];}
        }
    }
}
END { print "'${prefix}'\t"count"\t"len >> "similar/'${prefix}'_similar_count.txt" }' > similar/${prefix}.contig_inheritance_v0.8.postprocess.similar