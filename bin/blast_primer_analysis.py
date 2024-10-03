#!/usr/bin/env python

import pandas as pd
import os
import subprocess

def blast_sequence(primer_seq, index, reference_geno):
    # Output each sequence for BLAST
    sequence_id = f">{index}\n{primer_seq}\n"
    
    # Write temp sequence file
    with open("temp_sequence.fa", "w") as temp_seq_file:
        temp_seq_file.write(sequence_id)

    # Get the number of available CPU cores and subtract 1
    num_threads = max(os.cpu_count() - 1, 1)  # Ensure at least 1 thread

    # Run BLAST
    blast_command = [
        "blastn", 
        "-task", "blastn-short",
        "-query", "temp_sequence.fa",
        "-db", reference_geno,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-perc_identity", "100",
        "-qcov_hsp_perc", "100",
        "-num_threads", str(num_threads)  # Use available CPU cores minus one
    ]

    result = subprocess.run(blast_command, capture_output=True, text=True)

    return index, result.stdout.strip()

def main(reference_geno, output_file):
    # Read the filtered sequences from a file (e.g., "tmp_filtered")
    filtered_data = pd.read_csv("tmp_filtered", sep="\t")

    with open(output_file, "w") as output:
        # Write header to output file
        header = "index\tproduct_size\ttype\tstart\tend\tvariation_number\t3'diffall\tlength\ttm\tgccontent\tany\t3'\tend_stability\thairpin\tprimer_seq\treversecomplement\tpenalty\tcompl_any\tcompl_end\tscore\tstandardized_name\tref/alt\tdesigned_primer\tstrand_orientation\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n"
        output.write(header)

        # Loop through filtered sequences and run BLAST
        for index, row in filtered_data.iterrows():
            primer_seq = row['primer_seq']

            # Run BLAST for the primer sequence
            index, blast_output = blast_sequence(primer_seq, index, reference_geno)

            # Process BLAST output
            if blast_output:
                for line in blast_output.strip().splitlines():
                    # Split the BLAST output into columns
                    blast_columns = line.split("\t")
                    output_row = f"{index}\t{row['product_size']}\t{row['type']}\t{row['start']}\t{row['end']}\t{row['variation_number']}\t_\t{row['length']}\t{row['tm']}\t{row['gccontent']}\t{row['any']}\t{row['3\'']}\t{row['end_stability']}\t{row['hairpin']}\t{row['primer_seq']}\t{row['reversecomplement']}\t{row['penalty']}\t{row['compl_any']}\t{row['compl_end']}\t{row['score']}\t{row['standardized_name']}\t{row['ref/alt']}\t{row['designed_primer']}\t{row['strand_orientation']}\t" + "\t".join(blast_columns) + "\n"
                    output.write(output_row)
            else:
                # If no hits, append the original row with empty columns for BLAST results
                output_row = f"{index}\t{row['product_size']}\t{row['type']}\t{row['start']}\t{row['end']}\t{row['variation_number']}\t_\t{row['length']}\t{row['tm']}\t{row['gccontent']}\t{row['any']}\t{row['3\'']}\t{row['end_stability']}\t{row['hairpin']}\t{row['primer_seq']}\t{row['reversecomplement']}\t{row['penalty']}\t{row['compl_any']}\t{row['compl_end']}\t{row['score']}\t{row['standardized_name']}\t{row['ref/alt']}\t{row['designed_primer']}\t{row['strand_orientation']}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"
                output.write(output_row)

if __name__ == "__main__":
    import sys
    reference_geno = sys.argv[1]
    output_file = sys.argv[2]
    main(reference_geno, output_file)
