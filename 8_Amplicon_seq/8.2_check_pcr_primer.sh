for i in 4 22 30 36 51 52 61 77 96;do  ls callable_genome/individual/Bf-${i}/*homolog.list | while read line; do cat $line | sort | uniq | awk 'END{if(NR>0) {print "'${line}'""\t"NR}}';done ;done > homolog.list
cat homolog.list | sed 's/callable_genome\/individual\///g' | sed 's/_homolog.list//g' | cut -d \/ -f2 | sed 's/\./\t/g' > homolog.fmt.list
cat homolog.fmt.list | while read line; do {
    id=$(echo $line | cut -d " " -f1)
    contig=$(echo $line | cut -d " " -f2)
    homo=$(echo $line | cut -d " " -f3)
    fail_150=$(cat pcr_products/${id}/150/${id}.${contig}.len150.out | sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'BEGIN{fail=0}{sum=$5+$6+$7+$8;if(sum>2){fail+=1}}END{print fail}')
    fail_160=$(cat pcr_products/${id}/160/${id}.${contig}.len160.out | sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'BEGIN{fail=0}{sum=$5+$6+$7+$8;if(sum>2){fail+=1}}END{print fail}')
    fail_170=$(cat pcr_products/${id}/170/${id}.${contig}.len170.out | sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'BEGIN{fail=0}{sum=$5+$6+$7+$8;if(sum>2){fail+=1}}END{print fail}')
    fail_180=$(cat pcr_products/${id}/180/${id}.${contig}.len180.out | sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'BEGIN{fail=0}{sum=$5+$6+$7+$8;if(sum>2){fail+=1}}END{print fail}')
    fail_190=$(cat pcr_products/${id}/190/${id}.${contig}.len190.out | sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'BEGIN{fail=0}{sum=$5+$6+$7+$8;if(sum>2){fail+=1}}END{print fail}')
    fail_200=$(cat pcr_products/${id}/200/${id}.${contig}.len200.out | sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'BEGIN{fail=0}{sum=$5+$6+$7+$8;if(sum>2){fail+=1}}END{print fail}')
    home_150=$(cat pcr_products/${id}/150/${id}.${contig}.len150.out | sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'END{print NR}')
    home_160=$(cat pcr_products/${id}/160/${id}.${contig}.len160.out | sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'END{print NR}')
    home_170=$(cat pcr_products/${id}/170/${id}.${contig}.len170.out | sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'END{print NR}')
    home_180=$(cat pcr_products/${id}/180/${id}.${contig}.len180.out | sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'END{print NR}')
    home_190=$(cat pcr_products/${id}/190/${id}.${contig}.len190.out | sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'END{print NR}')
    home_200=$(cat pcr_products/${id}/200/${id}.${contig}.len200.out | sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'END{print NR}')
    awk -F "\t" 'BEGIN{OFS="\t";print "'${id}'","'${contig}'","'${homo}'","'${home_150}'","'${home_160}'","'${home_170}'","'${home_180}'","'${home_190}'","'${home_200}'","'${fail_150}'","'${fail_160}'","'${fail_170}'","'${fail_180}'","'${fail_190}'","'${fail_200}'"}'
} 
done | awk -F "\t" 'BEGIN{OFS="\t";print "ID\tcontig\thomo\t#150\t#160\t#170\t#180\t#190\t#200\tfail150\tfail160\tfail170\tfail180\tfail190\tfail200"}{OFS="\t";print $0}' > homolog_checked.fmt.list

cat homolog.fmt.list | while read line; do {
    id=$(echo $line | cut -d " " -f1)
    contig=$(echo $line | cut -d " " -f2)
    homo=$(echo $line | cut -d " " -f3)
    fail_150=$(cat pcr_products/${id}/150/${id}.${contig}.len150.pair1.out| sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'BEGIN{fail=0}{sum=$5+$6+$7+$8;if(sum>2){fail+=1}}END{print fail}')
    fail_160=$(cat pcr_products/${id}/160/${id}.${contig}.len160.pair1.out| sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'BEGIN{fail=0}{sum=$5+$6+$7+$8;if(sum>2){fail+=1}}END{print fail}')
    fail_170=$(cat pcr_products/${id}/170/${id}.${contig}.len170.pair1.out| sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'BEGIN{fail=0}{sum=$5+$6+$7+$8;if(sum>2){fail+=1}}END{print fail}')
    fail_180=$(cat pcr_products/${id}/180/${id}.${contig}.len180.pair1.out| sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'BEGIN{fail=0}{sum=$5+$6+$7+$8;if(sum>2){fail+=1}}END{print fail}')
    fail_190=$(cat pcr_products/${id}/190/${id}.${contig}.len190.pair1.out| sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'BEGIN{fail=0}{sum=$5+$6+$7+$8;if(sum>2){fail+=1}}END{print fail}')
    fail_200=$(cat pcr_products/${id}/200/${id}.${contig}.len200.pair1.out| sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'BEGIN{fail=0}{sum=$5+$6+$7+$8;if(sum>2){fail+=1}}END{print fail}')
    home_150=$(cat pcr_products/${id}/150/${id}.${contig}.len150.pair1.out| sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'END{print NR}')
    home_160=$(cat pcr_products/${id}/160/${id}.${contig}.len160.pair1.out| sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'END{print NR}')
    home_170=$(cat pcr_products/${id}/170/${id}.${contig}.len170.pair1.out| sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'END{print NR}')
    home_180=$(cat pcr_products/${id}/180/${id}.${contig}.len180.pair1.out| sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'END{print NR}')
    home_190=$(cat pcr_products/${id}/190/${id}.${contig}.len190.pair1.out| sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'END{print NR}')
    home_200=$(cat pcr_products/${id}/200/${id}.${contig}.len200.pair1.out| sed '1d' | awk -F "\t" '!array[$2]++' | awk -F "\t" 'END{print NR}')
    awk -F "\t" 'BEGIN{OFS="\t";print "'${id}'","'${contig}'","'${homo}'","'${home_150}'","'${home_160}'","'${home_170}'","'${home_180}'","'${home_190}'","'${home_200}'","'${fail_150}'","'${fail_160}'","'${fail_170}'","'${fail_180}'","'${fail_190}'","'${fail_200}'"}'
} 
done | awk -F "\t" 'BEGIN{OFS="\t";print "ID\tcontig\thomo\t#150\t#160\t#170\t#180\t#190\t#200\tfail150\tfail160\tfail170\tfail180\tfail190\tfail200"}{OFS="\t";print $0}' > homolog_checked.fmt.pair1.list
