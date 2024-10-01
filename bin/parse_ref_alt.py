#!/usr/bin/env python3

# Import pandas and sys
import pandas as pd
import sys

# Function to compute the complement of a nucleotide
def complement(base):
    complement_dict = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    return complement_dict.get(base.upper(), base)  # Return the complement or the base itself if not found

# Function to determine if the primer is Ref or Alt
def label_allele(primer_seq, ref_allele, alt_allele):
    last_base = primer_seq[-1]  # Get the last base of the primer sequence
    complement_base = complement(last_base)  # Get the complement of the last base
    
    # First check the last base
    if last_base == ref_allele:
        return "Ref"
    elif last_base == alt_allele:
        return "Alt"
    
    # If no match, check the complement
    if complement_base == ref_allele:
        return "Ref"
    elif complement_base == alt_allele:
        return "Alt"
    
    # If no matches, return Unknown
    return "Unknown"

def format_pos(pos):
    return f"{int(pos):09d}"  # Format pos to always be a 9-digit number

def main(primers_file, variants_file, standardized_prefix, output_file):
    # Load the primer and variant data from files
    primers_df = pd.read_csv(primers_file, sep="\t")
    variants_df = pd.read_csv(variants_file, sep="\t")

    # Replace all "_" with "-" in the variant IDs to match primer format
    variants_df['id'] = variants_df['id'].str.replace('_', '-', regex=False)

    # Create the standardized_name column in variants_df
    variants_df['standardized_name'] = standardized_prefix + '_' + variants_df['chr'] + variants_df['pos'].apply(format_pos)

    # Print for debugging
    print("Primers DataFrame:")
    print(primers_df.head())
    print("\nVariants DataFrame:")
    print(variants_df.head())

    # VIC and FAM sequences
    FAM = "GAAGGTGACCAAGTTCATGCT"
    VIC = "GAAGGTCGGAGTCAACGGATT"

    # Initialize variables to store the strand orientations
    ref_strand_orientation = ''
    alt_strand_orientation = ''

    # Loop through the primer data and cross-reference with variant data
    for index, row in primers_df.iterrows():
        # Extract primer ID (assuming format like S1A-3618173-left-0-Allele-A)
        primer_id = "-".join(row['index'].split("-")[:2])  # Get the form "S1A-3618173"
        
        # Use exact case-sensitive match to find a matching variant row
        variant_row = variants_df[variants_df['id'] == primer_id]
        
        # Initialize variables
        label = "Unknown"
        designed_primer = ""
        strand_orientation = ''
        
        # Check if there is a matching variant
        if not variant_row.empty:
            ref_allele = variant_row['ref'].values[0]
            alt_allele = variant_row['alt'].values[0]

            # Check if the primer is labeled as 'common'
            if 'common' in row['index'].lower():  # Case insensitive check
                label = 'CR'  # Change to "CR"
                designed_primer = row['primer_seq']  # Keep the common primer as is
                strand_orientation = '-'  # Set to default for common initially
            else:
                # Determine if the primer is Ref or Alt
                label = label_allele(row['primer_seq'], ref_allele, alt_allele)

                # Determine strand orientation and designed primer
                if label == "Ref":
                    strand_orientation = '+' if row['primer_seq'][-1] == ref_allele else '-'
                    designed_primer = FAM + row['primer_seq']  # Add FAM for Ref
                    ref_strand_orientation = strand_orientation  # Store Ref strand orientation
                elif label == "Alt":
                    strand_orientation = '+' if row['primer_seq'][-1] == alt_allele else '-'
                    designed_primer = VIC + row['primer_seq']  # Add VIC for Alt
                    alt_strand_orientation = strand_orientation  # Store Alt strand orientation

            # Add standardized_name from variant_row to primers_df
            primers_df.loc[index, 'standardized_name'] = variant_row['standardized_name'].values[0]

        # Modify the new labeled index column
        primers_df.loc[index, 'Ref/Alt'] = label
        primers_df.loc[index, 'Designed Primer'] = designed_primer
        primers_df.loc[index, 'Strand Orientation'] = strand_orientation

        # Update common primers based on Ref and Alt orientations
        if 'common' in row['index'].lower():  # Check for common primers
            # Determine common strand orientation based on Ref and Alt orientations
            if ref_strand_orientation == '-' or alt_strand_orientation == '-':
                primers_df.loc[index, 'Strand Orientation'] = '-'  # Match the strand orientation to Ref or Alt
            else:
                primers_df.loc[index, 'Strand Orientation'] = '+'  # Set to positive if both Ref and Alt are positive

    # Print for debugging after assigning standardized names
    print("\nDataFrame with Standardized Names:")
    print(primers_df[['index', 'standardized_name']].head())

    # Drop Labeled Index if it exists
    if 'Labeled Index' in primers_df.columns:
        primers_df.drop(columns=['Labeled Index'], inplace=True)

    # Clean column names: replace spaces with underscores and convert to lowercase
    primers_df.columns = primers_df.columns.str.replace(' ', '_').str.lower()

    # Save the modified primer data to a new .tsv file
    primers_df.to_csv(output_file, sep="\t", index=False)
    print(f"Labeled primer data saved to {output_file}")

if __name__ == "__main__":
    # Ensure that the correct number of arguments are provided
    if len(sys.argv) != 5:
        print("Usage: python label_primers.py <primers_file> <variants_file> <standardized_prefix> <output_file>")
        sys.exit(1)

    # Get the command line arguments
    primers_file = sys.argv[1]
    variants_file = sys.argv[2]
    standardized_prefix = sys.argv[3]  # Prefix for standardized_name
    output_file = sys.argv[4]  # Output filename

    # Run the main function with the provided file paths
    main(primers_file, variants_file, standardized_prefix, output_file)
