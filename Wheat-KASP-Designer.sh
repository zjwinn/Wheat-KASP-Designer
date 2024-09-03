#!/bin/bash

# Function to display script usage
usage() {
    echo
    echo "###############################################"
    echo "#                                             #"
    echo "#               Help Information              #"
    echo "#                                             #"
    echo "###############################################"
    echo
    echo "Usage:"
    echo -e "\t$0 [OPTIONS] ARGUMENT"
    echo
    echo "Description:"
    echo -e "\tThis script will take a known reference genome (.fa) and variant calling format (VCF) file"
    echo -e "\tand report back subgenome specific KASP markers for wheat (Triticum aestivum L.)"
    echo -e "\twith some metric for quality assessment."
    echo 
    echo "Options:"
    echo -e "\t-h, --help              Display this help and exit"
    echo -e "\t-v, --verbose           Display text feedback (default option is false)"
    echo 
    echo "Arguments:"
    echo -e "\t-i, --input-file        realpath to a properly formatted, gziped VCF file (.vcf.gz)"
    echo -e "\t-o, --output-file       name of the output file (tab-delimited)"
    echo -e "\t-r, --reference-geno    realpath to a reference genome file (.fa)"
    echo -e "\t-s, --snps              a .txt file with a vector of SNP names to subset from the VCF"
    echo
    echo "Examples:"
    echo -e "\t$(realpath $0) -i 'input_file.vcf' -o 'output.txt' -r 'reference_genome.fa'"
    exit 1
}

# Set verbose to false automatically
verbose=false

# Check if there are no argments
if [ $# -eq 0 ]; then
    echo "Error: No inputs provided" >&2
    usage
fi

# Parse command-line options
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input-file)
            input_file="$2"
            shift
            shift
            ;;
        -o|--output-file)
            output_file="$2"
            shift
            shift
            ;;            
        -r|--reference-geno)
            reference_geno="$2"
            shift
            shift
            ;;
        -s|--snps)
            snp_list="$2"
            shift 
            shift 
            ;;
        -h|--help)
            usage
            ;;
        -v|--verbose)
            verbose=true
            shift
            shift            
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check if required options are provided
if [ -z "$input_file" ] || [ -z "$output_file" ] || [ -z "$reference_geno" ]; then
    echo "Error: Input file, output file, and reference genome file are required." >&2
    usage
fi

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found." >&2
    usage
fi

# Check if reference genome
if [ ! -f "$reference_geno" ]; then
    echo "Error: Input file '$reference_geno' not found." >&2
    usage
fi

# Get realpath of snp file
if [ -n "$snp_list" ]; then
    # Convert to real path
    real_path_snp_list=$(realpath "$snp_list")
fi

if [ "$verbose" = true ]; then
    # Print header
    echo
    echo "###############################################"
    echo "#                                             #"
    echo "#          Wheat KASP Designer v1.0           #"
    echo "#                                             #"
    echo "###############################################"
    echo
    echo "Written by: Zachary J. Winn PhD"
    echo "Contact information:"
    echo -e "\tGovernment Email: zachary.winn@usda.gov"
    echo -e "\tPersonal Email: zwinn@outlook.com"
    echo
    echo "###############################################"
    echo "# WARNING: This program is not under warranty #"
    echo "#          Use at your own discretion!        #"
    echo "###############################################"
fi

# Create a temporary directory
tmp_dir=$(mktemp -d -t Wheat-KASP-Designer-XXXXXX)

# # Define a cleanup function
# cleanup() {
#     # Display message
#     if [ "$verbose" = true ]; then
#         echo
#         echo "### Cleaning up temporary directory: $tmp_dir"
#     fi
    
#     # Remove directory
#     rm -rf "$tmp_dir"

#     # Display message
#     if [ "$verbose" = true ]; then
#         echo
#         echo "### Temporary directory cleaned up. Exiting script."
#     fi
# }

# # Ensure the cleanup function is called on script exit
# trap cleanup EXIT

# Get working directory
working_directory=$(pwd)

# Pull the script location
script_dir=$(realpath $(dirname $0))

# Change to temp directory
cd $tmp_dir

# Debug info
echo
echo "### Temporary directory: $tmp_dir"
echo "### SNP list: $real_path_snp_list"
echo "### Current WD: $(pwd)"

# First filter .vcf
if [ -z "$snp_list" ]; then
    # Display message
    if [ "$verbose" = true ]; then
        echo
        echo "### No SNP list detected. Subsetting to only biallelic SNP."
    fi

    # Run vcftools to subset SNP
    vcftools --gzvcf "$input_file" \
        --min-alleles 2 \
        --max-alleles 2 \
        --temp "$temp_dir" \
        --recode \
        --stdout | bgzip > "$tmp_dir/filt.vcf.gz"
else
    # Display message
    if [ "$verbose" = true ]; then
        echo
        echo "### SNP list detected. Subsetting to SNPs indicated."
    fi

    # Run vcftools to subset SNP
    vcftools --gzvcf "$input_file" \
        --snps "$real_path_snp_list" \
        --min-alleles 2 \
        --max-alleles 2 \
        --temp "$temp_dir" \
        --recode \
        --stdout | bgzip > "$tmp_dir/filt.vcf.gz"
fi

# Now pull the chr, pos, id, ref, alt
echo -e "chr\tpos\tid\tref\talt" > snp_seq_pull_input.txt
zcat filt.vcf.gz | grep -v '^#' | cut -f1-5 >> snp_seq_pull_input.txt

# Do check
check1=$(zcat "$reference_geno" | head -n 1 | grep "Chr")

# Check if chromosome names have Chr in them 
if [ -n "$check1" ]; then
    # Process the file with awk, skip the first line, and write to a temporary file
    awk -F'\t' 'BEGIN {OFS="\t"} {if (NR > 1) $1="Chr"$1; print}' snp_seq_pull_input.txt > temp_file

    # Move the temporary file back to the original file
    mv temp_file snp_seq_pull_input.txt
fi

# Now get the sequences using snp_sequence_puller.sh
bash "$script_dir/SNP-Sequence-Puller/snp_sequence_puller_auto.sh" \
    -i ./snp_seq_pull_input.txt \
    -o ./snp_seq_pull_output.txt \
    -r "$reference_geno"

# Remove any error lines from the output of the above script
awk -F'\t' '$6 !~ /Error/' snp_seq_pull_output.txt > temp_file
mv temp_file snp_seq_pull_output.txt

# Change back to working directory
cd "$working_directory"

# Exit without error
exit 0
