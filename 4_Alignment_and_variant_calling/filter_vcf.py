#!/usr/bin/env python3
import pysam  # type: ignore
import argparse

def check_alt_in_bam(bam, chrom, pos, alt_allele):
    """
    Check if any read in BAM file at chromosome:position (VCF coordinates, 1-indexed) 
    contains the alt_allele. Uses pileup which expects 0-indexed positions, hence pos-1.
    """
    count = 0
    for pileupcolumn in bam.pileup(chrom, pos - 1, pos, truncate=True):
        # pileupcolumn.pos is 0-indexed
        if pileupcolumn.pos == pos - 1:
            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del or pileupread.is_refskip:
                    continue
                query_pos = pileupread.query_position
                if query_pos is None:
                    continue
                base = pileupread.alignment.query_sequence[query_pos]
                if base.upper() == alt_allele.upper():
                    count += 1
    return count >= 2

def main():
    parser = argparse.ArgumentParser(
        description="Filter VCF sites by checking if alt alleles are supported in both BAMs"
    )
    parser.add_argument("--in", dest="vcf_file", required=True, 
                       help="Input VCF path (must be bgzipped and indexed)")
    parser.add_argument("--bam", dest="bam_files", action="append", required=True,
                       help="BAM file path (provide twice)")
    parser.add_argument("--out", dest="out_vcf", required=True, 
                       help="Output VCF path")
    args = parser.parse_args()

    if len(args.bam_files) != 2:
        parser.error("Exactly two --bam files must be provided.")

    try:
        bam1 = pysam.AlignmentFile(args.bam_files[0], "rb")
        bam2 = pysam.AlignmentFile(args.bam_files[1], "rb")
    except Exception as e:
        parser.error(f"Failed to open BAM file: {e}")

    try:
        vcf_in = pysam.VariantFile(args.vcf_file)
    except Exception as e:
        parser.error(f"Failed to open VCF file: {e}")

    try:
        vcf_out = pysam.VariantFile(args.out_vcf, "w", header=vcf_in.header)
    except Exception as e:
        parser.error(f"Failed to create output VCF file: {e}")

    for record in vcf_in:
        if not record.alts:
            continue

        # Process only the first alt allele; can be extended to loop through all alts
        alt_allele = record.alts[0]

        found_bam1 = check_alt_in_bam(bam1, record.chrom, record.pos, alt_allele)
        found_bam2 = check_alt_in_bam(bam2, record.chrom, record.pos, alt_allele)

        # Write record to output if neither BAM supports the alt
        if not found_bam1 and not found_bam2:
            vcf_out.write(record)

    bam1.close()
    bam2.close()
    vcf_in.close()
    vcf_out.close()

if __name__ == "__main__":
    main()