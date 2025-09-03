import sys
import csv

def main():
    prefix = sys.argv[1]
    
    # Part 1: Generate similar_relation_between_parents.tsv
    ids = set()
    same = {}
    same_P = {}
    input_file = f"/public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/cosin_similarity/postprocess/{prefix}.contig_inheritance_v0.8.postprocess.out"
    output_tsv = f"tmp/{prefix}.similar_relation_between_parents.tsv"
    filtered = []        # Second filter condition
    with open("/public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/similar_conitgs/similar_relation_between_parents.tsv") as fref, open(input_file) as fin:        
        for line in csv.reader(fref, delimiter='\t'):
            if line[0] in same:
                same[line[0]].append(line[1])
            else:
                same[line[0]] = [line[1]]
            if line[1] in same_P:
                same_P[line[1]].append(line[0])
            else:
                same_P[line[1]] = [line[0]]
        for row in csv.reader(fin, delimiter='\t'):
            if not (row[0].startswith('non') and float(row[15]) > 0.5) and row[18].startswith('PASS'):
                filtered.append((row[0], row[16]))
                ids.add(row[0])
        for row in csv.reader(fin, delimiter='\t'):
            if row[0].startswith('non') and float(row[15]) > 0.5 and row[18].startswith('PASS'):
                if row[0] in same and row[0].endwith('M'):
                    for i in same[row[0]]:
                        if i in ids and row[0] not in ids:
                            ids.add(row[0])
                if row[0] in same_P and row[0].endwith('P'):
                    for i in same_P[row[0]]:
                        if i in ids and row[0] not in ids:
                            ids.add(row[0])
    with open(output_tsv, 'w') as fout, open("/public5/home/xuejing/Amphioxus/new_align/reassembly/read_align/similar_conitgs/similar_relation_between_parents.tsv") as fref:
        tsv_writer = csv.writer(fout, delimiter='\t')     
        for row in csv.reader(fref, delimiter='\t'):
            if (row[0] in ids) and (row[1] in ids):
                tsv_writer.writerow(row)

    # Part 2: Generate contig_inheritance_v0.8.postprocess.list
    output_list = f"bed/{prefix}.contig_inheritance_v0.8.postprocess.remove_similarity.bed"
    a = {}
    b = {}
    
    # Read temporary file to establish mapping relationship
    with open(output_tsv) as fin:
        for row in csv.reader(fin, delimiter='\t'):
            a[row[0]] = row[1]
            b[row[1]] = row[0]
    
    # Process final output
    with open(input_file) as fin, open(output_list, 'w') as fout:
        tsv_writer = csv.writer(fout, delimiter='\t')
        
        for row in csv.reader(fin, delimiter='\t'):
            if not (row[0].startswith('non') and float(row[15]) > 0.5) and row[18].startswith('PASS'):
                contig_id = row[0]
                parent = row[16]
                if (contig_id not in a) and (contig_id not in b):
                    tsv_writer.writerow([contig_id, 100, int(parent)-100])

if __name__ == '__main__':
    main()
