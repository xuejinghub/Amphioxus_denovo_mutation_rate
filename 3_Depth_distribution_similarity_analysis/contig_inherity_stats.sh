#!/bin/bash
#SBATCH -w node80
#SBATCH -c 1
#SBATCH --time=48:00:00                     
#SBATCH --mem=80G
#SBATCH --job-name=inheritance
#SBATCH --output=log/inheritance_%a.out


for i in 7;do cat /public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/cosin_similarity/postprocess/Bf-${i}.contig_inheritance_v0.8.postprocess.similar |  grep "PASS" | cut -f1,17,20 |  awk -F "\t" '{if($1~/M$/) {
    if($0~/^non/) {
        m_non+=1;m_non_len+=$2;if($3~"parental_Similar"){
            m_sim+=1
            } else{
            m_dissim+=1
        }
    } else if ($0~/^p/) {
        m_pri+=1;m_pri_len+=$2;if($3~"parental_Similar"){
        m_sim+=1
        } else{
            m_dissim+=1
        }
    } else if($0~/^s/) {
        m_sec+=1;m_sec_len+=$2;if($3~"parental_Similar"){
            m_sim+=1
        } else{
            m_dissim+=1
        }
    }
} else if($1~/P$/) {
    if($0~/^non/) {
        p_non+=1;p_non_len+=$2;if($3~"parental_Similar"){
            p_sim+=1
        } else{
            p_dissim+=1
        }
    } else if ($0~/^p/) {
        p_pri+=1;p_pri_len+=$2;if($3~"parental_Similar"){
            p_sim+=1
        } else{
            p_dissim+=1
        }
    } else if($0~/^s/) {
        p_sec+=1;p_sec_len+=$2;if($3~"parental_Similar"){
            p_sim+=1
        } else{
            p_dissim+=1
        }
    }
}}END{OFS="\t";print "Bf-'${i}'",m_non,m_non_len,m_pri,m_pri_len,m_sec,m_sec_len,p_non,p_non_len,p_pri,p_pri_len,p_sec,p_sec_len,m_non_len+m_pri_len+m_sec_len,p_non_len+p_pri_len+p_sec_len,(m_non_len+m_pri_len+m_sec_len)/(p_non_len+p_pri_len+p_sec_len),m_sim,m_dissim,p_sim,p_dissim}';done > test.bias7.tsv