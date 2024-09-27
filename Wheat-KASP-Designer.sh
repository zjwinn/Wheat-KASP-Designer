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
    echo -e "\tand report back subgenome-specific KASP markers for wheat (Triticum aestivum L.)"
    echo -e "\twith some metric for quality assessment."
    echo 
    echo "Options:"
    echo -e "\t-h, --help              Display this help and exit"
    echo -e "\t-v, --verbose           Display text feedback (default option is false)"
    echo -e "\t-a, --keep-anyway       indicates to keep markers even when failing filters" 
    echo -e "\t-k, --kasp              indicates to design KASP assays"
    echo -e "\t-c, --caps              indicates to design CAPS assays"     
    echo 
    echo "Arguments:"
    echo -e "\t-i, --input-file        realpath to a properly formatted, gzipped VCF file (.vcf.gz)"
    echo -e "\t-r, --reference-geno    realpath to a reference genome file (.fa)"
    echo -e "\t-s, --snps              a .txt file with a vector of SNP names to subset from the VCF"
    echo -e "\t-m, --max-temp          a maximum temperature indicated in Celsius (default 64)"
    echo -e "\t-p, --max-price         a maximum price in USD (default 200)"
    echo -e "\t-x, --max-size          a maximum size of the primer in base pairs (default 25)"
    echo -e "\t-j, --product-size      a maximum amplicon product size in base pairs (default 150)"  
    echo
    echo "Examples:"
    echo -e "\t$(realpath $0) -i 'input_file.vcf.gz' -o 'output.txt' -r 'reference_genome.fa' -s 'snps.txt' -k -m 63 -p 200 -x 25 -v"
    exit 1
}

# Set default options
verbose=false
debug=false
design_kasp=false
design_caps=false
max_temp=64
max_price=200
max_size=25
product_size_threshold=150
pick_anyway=0

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
        -m|--max-temp)
            max_temp="$2"
            shift 2
            ;;
        -p|--max-price)
            max_price="$2"
            shift 2
            ;;
        -x|--max-size)
            max_size="$2"
            shift 2
            ;;
        -j|--product-size)
            product_size_threshold="$2"
            shift 2
            ;;            
        -a|--keep-anyway)
            keep_anyway=1
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
    echo
    echo "######################################################"
    echo "### WARNING: DEBUG ACTIVE AND MAY MESS UP OUTPUTS! ###"
    echo "######################################################"
fi

# Check if required options are provided
if [ -z "$input_file" ] || [ -z "$reference_geno" ]; then
    echo "*** Error: Input file and reference genome file are required." >&2
    usage
fi

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "*** Error: Input file '$input_file' not found." >&2
    usage
fi

# Check if reference genome
if [ ! -f "$reference_geno" ]; then
    echo "*** Error: Input file '$reference_geno' not found." >&2
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
        # Remove directory
        rm -rf "$tmp_dir"

        # Display message
        if [ "$verbose" = true ]; then
            echo
            echo "#######################################################"
            echo "### Temporary directory cleaned up. Exiting script. ###"
            echo "#######################################################"
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
        echo "##################################################################################"
        echo "### BLAST database for the reference genome does not exist. Making database... ###"
        echo "##################################################################################"
    fi

    # Make the database
    makeblastdb -in "$reference_geno" -dbtype nucl -out "$reference_geno"
fi

if [ ! -f "$reference_geno_location/${reference_geno_file}.fai" ]; then
    # Display message
    if [ "$verbose" = true ]; then
        echo
        echo "#######################################################"
        echo "### Genome index file does not exist. Making index... #"
        echo "#######################################################"
    fi

    # Make the database
    samtools faidx "$reference_geno"
fi

# First filter .vcf
if [ -z "$snp_list" ]; then
    # Display message
    if [ "$verbose" = true ]; then
        echo
        echo "###############################################################"
        echo "### No SNP list detected. Subsetting to only biallelic SNP... #"
        echo "###############################################################"        
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
        echo "########################################################"          
        echo "### SNP list detected. Subsetting to SNPs indicated... #"
        echo "########################################################"
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

# Check verbose
if [ "$verbose" = true ]; then
    # Print
    echo
    echo "##################################"
    echo "# Generating polymarker input... #"
    echo "##################################"
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
python3 "$script_dir/SNP_Primer_Pipeline2/bin/parse_polymarker_input.py" snp_seq_pull_output.csv > parse_polymarker_input.py.log

# Check verbose
if [ "$verbose" = true ]; then
    # Print
    echo "###############################"
    echo "# BLASTing input sequences... #"
    echo "###############################"
fi

# Now blast
blastn -task blastn -db $reference_geno -query for_blast.fa -outfmt "6 std qseq sseq slen" -num_threads 8 -out blast_out.txt

# Check verbose
if [ "$verbose" = true ]; then
    # Print
    echo
    echo "###########################"
    echo "# Parsing BLAST output... #"
    echo "###########################"
fi

# Parse the blast output file and output the homelog contigs and flanking ranges
python3 "$script_dir/SNP_Primer_Pipeline2/bin/getflanking.py" snp_seq_pull_output.csv blast_out.txt temp_range.txt 3 > getflanking.py.log

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

# Check to make KASP 
if [ $design_kasp = true ]; then
    # Check verbose
    if [ "$verbose" = true ]; then
        # Print
        echo
        echo "###############################"
        echo "# Designing potential KASP... #"
        echo "###############################"
    fi    
    
    # Run design
    python3 "$script_dir/SNP_Primer_Pipeline2/bin/getkasp3.py" "$max_temp" "$max_size" "$pick_anyway" > getkasp3.py.log
fi

# Check to make CAPS
if [ $design_caps = true ]; then
    # Check verbose
    if [ "$verbose" = true ]; then
        # Print
        echo
        echo "###############################"
        echo "# Designing potential CAPS... #"
        echo "###############################"
    fi    
    
    # Run design
    python3 "$script_dir/SNP_Primer_Pipeline2/bin/getCAPS.py" "$max_price" "$max_temp" "$max_size" "$pick_anyway" > getCAPS.py.log
fi

# Define the expected number of columns
expected_columns=20

# Function to filter files based on column count and concatenate them
concat_filtered_files() {
    local input_dir=$1
    local output_file=$2
    local file_pattern=$3

    if [ "$debug" = true ]; then
        echo
        echo "#############################"
        echo "# File parsing script debug #"
        echo "#############################"
        echo 
        echo "### Processing files in directory: $input_dir"
        echo "### Looking for files matching pattern: $file_pattern"
        echo
    fi

    # Get the header from the first file and store it in the output file
    first_file=$(find $input_dir -name "$file_pattern*" | head -n 1)
    
    if [ "$debug" = true ]; then
        echo "### Using header from file: $first_file"
        echo
    fi

    head -n 1 "$first_file" > $output_file

    # Loop through files in the directory
    find $input_dir -name "$file_pattern*" | while read file; do
        if [ "$debug" = true ]; then
            echo "#### Processing file: $file"
            echo
        fi

        # Skip the header in each file and filter rows based on the expected number of columns
        tail -n +2 "$file" | while read line; do
            # Print the line and count the number of fields (columns)
            num_fields=$(echo "$line" | awk -F'\t' '{print NF}')

            if [ "$debug" = true ]; then
                echo "### Line: $line"
                echo "### Number of columns: $num_fields"
            fi

            # If the line has the expected number of columns, append it to the output file
            if [ "$num_fields" -eq "$expected_columns" ]; then
                echo "$line" >> $output_file
                if [ "$debug" = true ]; then echo; fi
            else
                if [ "$debug" = true ]; then
                    echo "### Skipping line due to incorrect number of columns."
                    echo
                fi
            fi
        done
    done
}

# Check verbose
if [ "$verbose" = true ]; then
    # Print
    echo
    echo "#######################################################"
    echo "# Concating output and performing BLAST per primer... #"
    echo "#######################################################"
fi 

# Check if there is KASP output and concatenate if there is (silently)
if [ $(find ./KASP_output/ -name "selected_KASP_primers*" 2>/dev/null | wc -l) -gt 0 ]; then
    # Send debug messages
    if [ "$debug" = true ]; then
        echo "### Concating KASP files..."
        echo
    fi

    # Step 1: Concat
    concat_filtered_files "./KASP_output" "Potential_KASP_primers.tsv" "selected_KASP_primers"
    
    # Send debug messages
    if [ "$debug" = true ]; then
        echo "### Filtering KASP files..."
        echo
    fi

    # Step 2: Filter based on product size
    awk -v threshold="$product_size_threshold" '$2 > threshold' "Potential_KASP_primers.tsv" > "tmp_filtered"
    awk -v threshold="$product_size_threshold" -F"\t" '{if ($2 + 0 <= threshold) print $0;}' "Potential_KASP_primers.tsv" > "tmp_filtered"


    # Send debug messages
    if [ "$debug" = true ]; then
        echo "### Creating template for BLAST of KASP files..."
        echo
    fi

    # Step 3: Prepare the output file with headers
    echo -e "index\tproduct_size\ttype\tstart\tend\tvariation\t3'diffall\tlength\tTm\tGCcontent\tany\t3'\tend_stability\thairpin\tprimer_seq\tReverseComplement\tpenalty\tcompl_any\tcompl_end\tscore\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > "tmp_blast_output.tsv"
    
    # Skip first line
    first_line=true

    # Step 4: Loop through filtered sequences and run BLAST
    while IFS=$'\t' read -r index product_size type start end variation _ length tm gc_content any three_end end_stability hairpin primer_seq reverse_complement penalty compl_any compl_end score; do
        # Skip first line
        if [ "$first_line" = true ]; then
            first_line=false
            continue  # Skip the first line
        fi       
        
        # Output each sequence for BLAST
        echo ">${index}" > temp_sequence.fa
        echo "$primer_seq" >> temp_sequence.fa
            
        # Send debug messages
        if [ "$debug" = true ]; then
            echo "### BLASTING  KASP files line..."
            echo
        fi

        # Run BLAST (you can customize the parameters as needed)
        blastn -task blastn-short \
            -query temp_sequence.fa \
            -db $reference_geno \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -perc_identity 100 \
            -qcov_hsp_perc 100  > temp_blast.txt
        
        # Check if BLAST returned hits
        if [ -s temp_blast.txt ]; then
            # Append each BLAST hit to the original row
            while IFS=$'\t' read -r qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore; do
                # Echo line 
                echo -e "${index}\t${product_size}\t${type}\t${start}\t${end}\t${variation}\t_\t${length}\t${tm}\t${gc_content}\t${any}\t${three_end}\t${end_stability}\t${hairpin}\t${primer_seq}\t${reverse_complement}\t${penalty}\t${compl_any}\t${compl_end}\t${score}\t$sseqid\t$pident\t$length\t$mismatch\t$gapopen\t$qstart\t$qend\t$sstart\t$send\t$evalue\t$bitscore" >> "tmp_blast_output.tsv"
            done < temp_blast.txt
        else
            # If no hits, append the original row with empty columns for BLAST results
            echo -e "${index}\t${product_size}\t${type}\t${start}\t${end}\t${variation}\t_\t${length}\t${tm}\t${gc_content}\t${any}\t${three_end}\t${end_stability}\t${hairpin}\t${primer_seq}\t${reverse_complement}\t${penalty}\t${compl_any}\t${compl_end}\t${score}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> "tmp_blast_output.tsv"
        fi
    done < "tmp_filtered"

    # Step 5: Move the BLAST-processed data to the final output
    mv "tmp_blast_output.tsv" "$working_directory/Potential_KASP_primers.tsv"
fi


# Check if there is CAPS output and concatenate if there is (silently)
if [ $(find ./CAPS_output/ -name "selected_CAPS_primers*" 2>/dev/null | wc -l) -gt 0 ]; then
    concat_filtered_files "./CAPS_output" "$working_directory/Potential_CAPS_primers.tsv" "selected_CAPS_primers"
fi

# Change back to working directory
cd "$working_directory"

# Check verbose
if [ "$verbose" = true ]; then
    # Print
    echo
    echo "######################"
    echo "# Pipeline complete! #" 
    echo "######################"
    echo
    echo "##############################################################################"
    echo "# Check for 'Potential_KASP_primers.tsv' and/or 'Potential_CAPS_primers.tsv' #" 
    echo "# in current working directory.                                              #"
    echo "##############################################################################"
fi  

# Exit without error
exit 0
