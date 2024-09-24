#!/bin/bash

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
    #echo -e "\t-d, --debug             Enable debug mode (default option is false)"
    echo 
    echo "Arguments:"
    echo -e "\t-i, --input-file        realpath to a properly formatted, gzipped VCF file (.vcf.gz)"
    echo -e "\t-o, --output-file       name of the output file (tab-delimited)"
    echo -e "\t-r, --reference-geno    realpath to a reference genome file (.fa)"
    echo -e "\t-s, --snps              a .txt file with a vector of SNP names to subset from the VCF"
    echo -e "\t-k, --kasp              indicates to design KASP assays"
    echo -e "\t-c, --caps              indicates to design CAPS assays"    
    echo
    echo "Examples:"
    echo -e "\t$(realpath $0) -i 'input_file.vcf.gz' -o 'output.txt' -r 'reference_genome.fa'"
    exit 1
}

# Set verbose, debug, and assay design options to false by default
verbose=false
debug=false
design_kasp=false
design_caps=false

# Check if there are no arguments
if [ $# -eq 0 ]; then
    echo "Error: No inputs provided" >&2
    usage
fi

# Parse command-line options
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -i|--input-file)
            input_file=$(realpath "$2")
            shift 2
            ;;
        -o|--output-file)
            output_file=$(realpath "$2")
            shift 2
            ;;
        -r|--reference-geno)
            reference_geno=$(realpath "$2")
            shift 2
            ;;
        -s|--snps)
            snp_list=$(realpath "$2")
            shift 2
            ;;
        -k|--kasp)
            design_kasp=true
            shift
            ;;
        -c|--caps)
            design_caps=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        -d|--debug)
            debug=true
            shift
            ;;
        -v|--verbose)
            verbose=true
            shift
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Enable debug mode if debug is true
if [ "$debug" = true ]; then
    echo "######################################################"
    echo "### WARNING: DEBUG ACTIVE AND MAY MESS UP OUTPUTS! ###"
    echo "######################################################"
fi

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

# Check verbose
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

# Get working directory     
working_directory=$(pwd)    

# Pull the script location  
script_dir=$(realpath $(dirname $0))

# Create a temporary directory
if [ "$debug" = true ]; then
    tmp_dir="$working_directory/tmp"
    mkdir $tmp_dir
else
    tmp_dir=$(mktemp -d -t Wheat-KASP-Designer-XXXXXX)
fi   

# Make a cleanup function if debug is false
if [ "$debug" = false ]; then
    # Define a cleanup function
    cleanup() {
        # Display message
        if [ "$verbose" = true ]; then
            echo
            echo "### Cleaning up temporary directory: $tmp_dir"
        fi
        
        # Remove directory
        rm -rf "$tmp_dir"

        # Display message
        if [ "$verbose" = true ]; then
            echo
            echo "### Temporary directory cleaned up. Exiting script."
        fi
    }

    # Ensure the cleanup function is called on script exit
    trap cleanup EXIT
fi

# Change to temp directory
cd $tmp_dir

# Debug info
if [ "$debug" = true ]; then
    echo
    echo "### Temporary directory: $tmp_dir"
    echo "### SNP list: $real_path_snp_list"
    echo "### Current WD: $(pwd)"
fi

# Get reference genome location and file name
reference_geno_location=$(dirname "$reference_geno")
reference_geno_file=$(basename "$reference_geno")

# Check if the database exists or not
if [ ! -f "$reference_geno_location/${reference_geno_file}.nal" ] && \
   [ ! -f "$reference_geno_location/${reference_geno_file}.ndb" ] && \
   [ ! -f "$reference_geno_location/${reference_geno_file}.njs" ]; then
    # Display message
    if [ "$verbose" = true ]; then
        echo
        echo "### BLAST database for the reference genome does not exist. Making database."
    fi

    # Make the database
    makeblastdb -in "$reference_geno" -dbtype nucl -out "$reference_geno"
fi

if [ ! -f "$reference_geno_location/${reference_geno_file}.fai" ]; then
    # Display message
    if [ "$verbose" = true ]; then
        echo
        echo "### Genome index file does not exist. Making index."
    fi

    # Make the database
    samtools faidx "$reference_geno"
fi

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
if [[ "$reference_geno" == *.gz ]]; then
    # Use zcat for .gz files
    check1=$(zcat "$reference_geno" | head -n 1 | grep "Chr")
else
    # Use cat for other files
    check1=$(cat "$reference_geno" | head -n 1 | grep "Chr")
fi

# Check if chromosome names have Chr in them 
if [ -n "$check1" ]; then
    # Process the file with awk, skip the first line, and write to a temporary file
    awk -F'\t' 'BEGIN {OFS="\t"} {if (NR > 1) $1="Chr"$1; print}' snp_seq_pull_input.txt > temp_file

    # Move the temporary file back to the original file
    mv temp_file snp_seq_pull_input.txt
fi

# Now get the sequences using snp_sequence_puller.sh
if [ "$verbose" = true ]; then
    # Run with verbose option
    bash "$script_dir/SNP-Sequence-Puller/snp_sequence_puller_auto.sh" \
        -i ./snp_seq_pull_input.txt \
        -o ./snp_seq_pull_output.txt \
        -r "$reference_geno" \
        -v
else
    # Run silent
    bash "$script_dir/SNP-Sequence-Puller/snp_sequence_puller_auto.sh" \
        -i ./snp_seq_pull_input.txt \
        -o ./snp_seq_pull_output.txt \
        -r "$reference_geno"
fi

# Remove any error lines from the output of the above script
awk -F'\t' '$6 !~ /Error/' snp_seq_pull_output.txt > temp_file
awk -F'\t' '$6 ~ /Error/' snp_seq_pull_output.txt > temp_file1
mv temp_file snp_seq_pull_output.txt
mv temp_file1 snp_seq_pull_output_errors.txt

# Check if BLAST-plus is installed
if ! command -v blastn &> /dev/null; then
    echo "### Error: BLAST+ is not installed."
    exit 1
fi

# Alter the file so that it looks correct for the python scripts
awk 'NR > 1 {print $3 "," $1 "," $6}' snp_seq_pull_output.txt > snp_seq_pull_output.csv

# Parse polymarker
python3 "$script_dir/bin/parse_polymarker_input.py" snp_seq_pull_output.csv > parse_polymarker_input.py.log

# Now blast
blastn -task blastn -db $reference_geno -query for_blast.fa -outfmt "6 std qseq sseq slen" -num_threads 8 -out blast_out.txt

# Parse the blast output file and output the homelog contigs and flanking ranges
python3 "$script_dir/bin/getflanking.py" snp_seq_pull_output.csv blast_out.txt temp_range.txt 3 > getflanking.py.log

# Replace - with \t in temp_range.txt
awk -F '\t' 'BEGIN {OFS = FS} {gsub(/-/, "\t", $3); print}' temp_range.txt > temp_range_mod.txt

# Get the size of the file
filesize=$(stat -c%s "temp_range.txt")  

# Check if the file size is 0
if [ "$filesize" -eq 0 ]; then  
    # Echo issue
    echo "All the SNPs are bad, possibly too many hits. Please check the stdout." | tee Potential_CAPS_primers.tsv Potential_KASP_primers.tsv > All_alignment_raw.fa
    
    # Exit the script with a non-zero status
    exit 1  
else
    # split file for each marker
    gawk -v OFS='\t' '{ print $2, $3, $4 > "temp_marker_" $1 ".bed" }' temp_range_mod.txt
fi

# For each marker in temp_marker
for i in temp_marker*; do 
    # Extract flanking sequence from database
    bedtools getfasta -fi $reference_geno -bed $i  -fo flanking_$i.fa
done

if [ $design_kasp = true ]; then
    cmd6 = script_path + "getkasp3.py " + max_Tm + " " + max_size + " " + pick_anyway # add blast option
    print("Step 6: Get KASP primers for each marker command:\n", cmd6)
    call(cmd6, shell=True)
fi

if [ $design_caps = true ]; then
    # Run design
    
fi

# Change back to working directory
cd "$working_directory"

# Exit without error
exit 0
