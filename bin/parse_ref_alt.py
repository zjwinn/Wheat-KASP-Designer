#!/usr/bin/env python3

# Import the necessary libraries
import pandas as pd
import sys

# Function to compute the complement of a nucleotide (used for matching reverse strands)
def complement(base):
    complement_dict = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    return complement_dict.get(base.upper(), base)  # Return the complement or the base itself if it's not A, T, C, or G

# Function to determine if a primer corresponds to the reference (Ref) or alternative (Alt) allele
def label_allele(primer_seq, ref_allele, alt_allele):
    last_base = primer_seq[-1]  # Get the last base of the primer sequence
    complement_base = complement(last_base)  # Get the complement of the last base

    # Check if the last base of the primer sequence matches the reference allele
    if last_base == ref_allele:
        return "Ref"
    # Check if the last base matches the alternative allele
    elif last_base == alt_allele:
        return "Alt"
    
    # If no match, check the complement of the last base for Ref or Alt
    if complement_base == ref_allele:
        return "Ref"
    elif complement_base == alt_allele:
        return "Alt"
    
    # If no match, return "Unknown"
    return "Unknown"

# Function to format the position column as a 9-digit number
def format_pos(pos):
    return f"{int(pos):09d}"  # Converts the position to an integer and pads it with zeros to be 9 digits long

def main(primers_file, variants_file, standardized_prefix, output_file):
    # Load the primer data and variant data from the provided files
    primers_df = pd.read_csv(primers_file, sep="\t")
    variants_df = pd.read_csv(variants_file, sep="\t")

    # Replace all underscores in the variant IDs with hyphens to match the format used in primers
    variants_df['id'] = variants_df['id'].str.replace('_', '-', regex=False)

    # Create a new standardized_name column by combining the prefix, chromosome, and formatted position
    variants_df['standardized_name'] = standardized_prefix + '_' + variants_df['chr'] + variants_df['pos'].apply(format_pos)

    # Define the fixed sequences for AL1 (FAM) and AL2 (VIC) primers
    FAM = "GAAGGTGACCAAGTTCATGCT"  # AL1
    VIC = "GAAGGTCGGAGTCAACGGATT"  # AL2

    # Initialize variables to store the strand orientations for reference and alternative alleles
    ref_strand_orientation = ''
    alt_strand_orientation = ''
    cr_counter = 1  # Initialize counter for common reverse primers (CR primers)

    last_allele_label = None  # Variable to track the last allele label (AL1 or AL2)
    last_orientation = '+'  # Variable to store the last assigned orientation for consistency

    # Loop through each row in the primer DataFrame
    for index, row in primers_df.iterrows():
        # Debug: Print the current primer row being processed
        print(f"Processing primer row {index}: {row}")

        # Extract the primer ID (assumes primer ID format like S1A-3618173)
        primer_id = "-".join(row['index'].split("-")[:2])  # Extract just the first two parts of the primer ID
        
        # Debug: Print the primer ID
        print(f"Primer ID: {primer_id}")
        
        # Find the corresponding variant in the variant DataFrame using the primer ID
        variant_row = variants_df[variants_df['id'] == primer_id]

        # Debug: Print the matched variant row
        if variant_row.empty:
            print(f"No matching variant found for primer ID: {primer_id}")
            continue  # Skip to the next iteration if no match found
        
        # Extract the reference and alternative alleles from the variant row
        ref_allele = variant_row['ref'].values[0]
        alt_allele = variant_row['alt'].values[0]

        # Initialize the standardized_name as a string from the variant_row
        standardized_name = variant_row['standardized_name'].values[0]

        # Initialize default values for primer labeling
        label = "Unknown"
        designed_primer = ""
        strand_orientation = ''
        common_reverse_label = ''  # For storing the common reverse label (CR1, CR2, etc.)
        
        # Check if the primer is a common reverse primer (contains 'common' in its name)
        if 'common' in row['index'].lower():
            label = 'CR'  # Label common primers as 'CR'
            designed_primer = row['primer_seq']  # Keep the common primer sequence as is
            # Set the strand orientation to the last assigned orientation for consistency
            primers_df.loc[index, 'Strand Orientation'] = last_orientation
            common_reverse_label = f"CR{cr_counter}"  # Add CR followed by the current counter value
            cr_counter += 1  # Increment the CR counter for the next common reverse primer
            
            # Update standardized_name for common reverse primers
            primers_df.loc[index, 'standardized_name'] = f"{standardized_name}-CR{cr_counter - 1}"
            
        else:
            # Reset the common reverse counter if the primer is labeled as Ref or Alt
            cr_counter = 1
            
            # Determine if the primer is for the Ref or Alt allele using the last base
            label = label_allele(row['primer_seq'], ref_allele, alt_allele)
            
            # Debug: Print the allele label
            print(f"Allele label for primer: {label}")
            
            # If the primer matches the reference allele
            if label == "Ref":
                strand_orientation = '+' if row['primer_seq'][-1] == ref_allele else '-'
                designed_primer = FAM + row['primer_seq']  # Add the FAM (AL1) sequence to the primer
                primers_df.loc[index, 'standardized_name'] = f"{standardized_name}-AL1"  # Update the standardized_name to label it as AL1
                
            # If the primer matches the alternative allele
            elif label == "Alt":
                strand_orientation = '+' if row['primer_seq'][-1] == alt_allele else '-'
                designed_primer = VIC + row['primer_seq']  # Add the VIC (AL2) sequence to the primer
                primers_df.loc[index, 'standardized_name'] = f"{standardized_name}-AL2"  # Update the standardized_name to label it as AL2

        # Set the Ref/Alt label, designed primer sequence, and strand orientation in the DataFrame
        primers_df.loc[index, 'Ref/Alt'] = label
        primers_df.loc[index, 'Designed Primer'] = designed_primer
        
        # Update the last assigned orientation after processing Ref or Alt primers
        if label == "Ref" or label == "Alt":
            last_allele_label = label  # Update the last allele label
            last_orientation = strand_orientation  # Store the last assigned orientation
        
        primers_df.loc[index, 'Strand Orientation'] = strand_orientation  # Assign the strand orientation for Ref/Alt primers

    # Clean up the DataFrame column names by converting spaces to underscores and lowercase
    primers_df.columns = primers_df.columns.str.replace(' ', '_').str.lower()

    # Save the updated primers DataFrame to the specified output file
    primers_df.to_csv(output_file, sep="\t", index=False)
    print(f"Labeled primer data saved to {output_file}")

# Entry point for the script
if __name__ == "__main__":
    # Ensure that the script is called with the correct number of arguments
    if len(sys.argv) != 5:
        print("Usage: python label_primers.py <primers_file> <variants_file> <standardized_prefix> <output_file>")
        sys.exit(1)

    # Read in command-line arguments
    primers_file = sys.argv[1]
    variants_file = sys.argv[2]
    standardized_prefix = sys.argv[3]
    output_file = sys.argv[4]

    # Call the main function to process the primers and variants
    main(primers_file, variants_file, standardized_prefix, output_file)
