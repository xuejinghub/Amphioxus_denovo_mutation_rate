import re
import csv
from collections import defaultdict

def process_contig_pairs(input_file, output_file):
    """Process raw result file and generate corrected output"""
    # Read raw data and build index
    records = []
    header_indices = {}
    with open(input_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)
        header_indices = {h: idx for idx, h in enumerate(headers)}
        records = list(reader)

    # Build grouping dictionary {bubbleNumber_source: {type: (row_index, metrics_dict)}}
    bubble_groups = defaultdict(lambda: defaultdict(dict))
    
    # First pass: build group index
    for row_idx, row in enumerate(records):
        contig_id = row[header_indices["ContigID"]]
        if m := re.match(r'^(primary|secondary)_bubble(\d+).*?([MP])$', contig_id):
            contig_type, bubble_num, source = m.groups()
            key = f"{bubble_num}_{source}"
            # Convert row data to metrics dictionary
            metrics = {header: row[idx] for header, idx in header_indices.items()}
            bubble_groups[key][contig_type] = (row_idx, metrics)

    # Second pass: resolve conflicts
    for group_key, contig_types in bubble_groups.items():
        if 'primary' in contig_types and 'secondary' in contig_types:
            # Get row indices and metrics data
            pri_idx, pri_met = contig_types['primary']
            sec_idx, sec_met = contig_types['secondary']
            
            # Convert to float for comparison
            pri_mp = float(pri_met['MP_hybrid'])
            sec_mp = float(sec_met['MP_hybrid'])
            
            # Determine target and conflict contigs
            if pri_mp < sec_mp:
                target_idx, conflict_idx = pri_idx, sec_idx
            else:
                target_idx, conflict_idx = sec_idx, pri_idx
            
            # Get original inheritance labels
            target_inheritance = records[target_idx][header_indices['inheritance']]
            conflict_inheritance = records[conflict_idx][header_indices['inheritance']]
            
            if target_inheritance == "PASS" and conflict_inheritance == "PASS":
                new_label = "conflict_resolved"
                records[conflict_idx][header_indices['inheritance']] = new_label
            # Handle dual non-PASS conflicts
            elif target_inheritance == "low_depth" and conflict_inheritance != "low_depth" and conflict_inheritance!= "PASS":
                new_label = "PASS_conflict_resolved"
                records[conflict_idx][header_indices['inheritance']] = new_label
            elif target_inheritance == "low_depth" and conflict_inheritance== "low_depth":
                continue
            elif target_inheritance!= "low_depth" and target_inheritance != "PASS" and conflict_inheritance== "low_depth":
                new_label = "PASS_conflict_resolved"
                records[target_idx][header_indices['inheritance']] = new_label
            elif target_inheritance != "PASS" and conflict_inheritance != "PASS":
                new_label = "PASS_conflict_resolved"
                records[conflict_idx][header_indices['inheritance']] = new_label

    # Write corrected file
    with open(output_file, 'w', newline='\n') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(headers)
        writer.writerows(records)
    print(f"Corrected results saved to: {output_file}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Resolve contig inheritance conflicts')
    parser.add_argument('-i', '--input', required=True, help='Input file path')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    
    args = parser.parse_args()
    process_contig_pairs(args.input, args.output)
