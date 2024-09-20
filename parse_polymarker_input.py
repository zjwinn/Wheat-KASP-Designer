#!/usr/bin/env python3
import sys

# Dictionary to map IUPAC codes to their corresponding single-letter codes
iupac = {
    "[A/G]": "R", "[G/A]": "R", "[C/T]": "Y", "[T/C]": "Y",
    "[G/C]": "S", "[C/G]": "S", "[A/T]": "W", "[T/A]": "W",
    "[G/T]": "K", "[T/G]": "K", "[A/C]": "M", "[C/A]": "M"
}

def main():
    # Check if the correct number of arguments is provided
    if len(sys.argv) < 2:
        print("Usage: python3 parse_polymarker_input.py <input_file> [--verbose]")
        sys.exit(1)

    # Get the input file and verbose flag from command-line arguments
    polymarker_input = sys.argv[1]
    verbose = "--verbose" in sys.argv

    # Open the output file for writing
    outfile = "for_blast.fa"
    out = open(outfile, "w")

    # Open the input file and skip the first line (header)
    with open(polymarker_input, "r") as infile:
        next(infile)  # Skip the header line

        # Process each line in the input file
        for line in infile:
            line = line.strip()  # Remove leading/trailing whitespace
            if not line:
                continue  # Skip empty lines

            # Print the original line if verbose mode is enabled
            if verbose:
                print(f"Original line: {line}")

            # Split the line into columns based on tab characters
            columns = line.split("\t")
            if verbose:
                print(f"Split columns: {columns}")

            # Ensure the line has exactly 6 columns
            if len(columns) != 6:
                if verbose:
                    print("Error: Line does not have exactly 6 columns")
                continue

            # Assign columns to variables
            chrom, pos, snpname, ref, alt, seq = columns

            # Replace underscores in SNP name with hyphens
            snpname = snpname.replace("_", "-")

            # Strip any extra spaces from the sequence
            seq = seq.strip()

            # Find the position of the SNP in the sequence
            pos = seq.find("[")

            # Print the content of variables if verbose mode is enabled
            if verbose:
                print(f"chrom: {chrom}, pos: {pos}, snpname: {snpname}, ref: {ref}, alt: {alt}, seq: {seq}")

            # Get the IUPAC code for the SNP
            snp = iupac.get(seq[pos:pos+5], 'Key not found')

            # Construct the new sequence with the IUPAC code
            seq2 = seq[:pos] + snp + seq[pos+5:]

            # Write the formatted output to the file
            out.write(f">{snpname}_{chrom}_{snp}\n{seq2}\n")

    # Close the output file
    out.close()

    return 0

if __name__ == '__main__':
    main()
