#!/usr/bin/env python3

import os
import glob
import sys

def concat_filtered_files(input_dir, output_file, file_pattern, expected_columns, debug=False):
    if debug:
        print("\n### File parsing script debug")
        print(f"### Processing files in directory: {input_dir}")
        print(f"### Looking for files matching pattern: {file_pattern}\n")

    first_file = glob.glob(os.path.join(input_dir, f"{file_pattern}*"))[0]
    
    if debug:
        print(f"### Using header from file: {first_file}\n")
    
    with open(first_file, 'r') as f:
        header = f.readline()
    
    with open(output_file, 'w') as f:
        f.write(header)  # Write the header to the output file
    
    for file in glob.glob(os.path.join(input_dir, f"{file_pattern}*")):
        if debug:
            print(f"#### Processing file: {file}\n")

        with open(file, 'r') as f:
            next(f)  # Skip the header
            for line in f:
                num_fields = len(line.strip().split('\t'))

                if debug:
                    print(f"### Line: {line.strip()}")
                    print(f"### Number of columns: {num_fields}")

                if num_fields == expected_columns:
                    with open(output_file, 'a') as out_f:
                        out_f.write(line)  # Append to the output file
                else:
                    if debug:
                        print("### Skipping line due to incorrect number of columns.\n")

if __name__ == "__main__":
    # Get command-line arguments
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    file_pattern = sys.argv[3]
    expected_columns = int(sys.argv[4])
    
    # Optional debug flag
    debug = '--debug' in sys.argv

    concat_filtered_files(input_dir, output_file, file_pattern, expected_columns, debug)
