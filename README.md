# Introduction
This repository contains a submodules of the popular [primer3](https://github.com/primer3-org/primer3) Git repository. In this repository, I (@zjwinn) will be writing a pipeline which may be used to take a properly formatted variant calling format (VCF) file and a reference genome (e.g., [Chinese Spring RefSeq V2.1](https://onlinelibrary.wiley.com/doi/10.1111/tpj.15289))output potential Kompetative allele specific primer ([KASP](https://www.biosearchtech.com/support/education/kasp-genotyping-reagents/how-does-kasp-work)) designs for biallelic single nucleotide polymoprhisms (SNPs). The purpose of this repository is to create an exicutable batch file (.sh) which can be called upon to design KASP automatically for mapping population or genome-wide marker panels. 

# Primer3 Information

Design PCR primers from DNA sequence. Widely used (190k Google hits for "primer3").
From mispriming libraries to sequence quality data to the generation of internal
oligos, primer3 does it. C&perl.

## Installing

```{bash}
`sudo apt-get install -y build-essential g++ cmake git-all`

`git clone https://github.com/primer3-org/primer3.git primer3`

`cd primer3/src`

`make`

`make test`
```

## Run Primer3

```{bash}
`./primer3_core ../example`
````

## Read the complete Primer3 manual

[Primer3 Manual](http://primer3.org/manual.html)

## SNP-Sequence-Puller Information

### Introduction
Single nucleotide polymorphisms (SNPs) are single nucleotide base-pair exchanges that occur genome wide in all living organisms. These SNPs can be used to design polymerase chain reaction (PCR) molecular markers. This program takes a reference genome, an alternate allele at a SNP site, and a specified flanking sequence length. The program then takes those inputs and outputs a sequence which is of some length with a centeralized SNP in the middle of the strand (I.E., ATTAG[C\T]GTACG). This output can then be taken and used in downstream process of molecular marker design. 

### Featured Scripts
There are two shell scripts written to acomplish the required task of pulling flanking sequences around SNPs and formatting them into an output:
1. [snp_sequence_puller.sh](https://github.com/zjwinn/SNP-Sequence-Puller/blob/main/snp_sequence_puller.sh)
2. [snp_sequence_puller_auto.sh](https://github.com/zjwinn/SNP-Sequence-Puller/blob/main/snp_sequence_puller_auto.sh)

The snp_sequence_puller.sh script takes a single provided SNP/position and returns a SNP with flanking sequence. The snp_sequence_puller_auto.sh script takes a tab delmited file of SNPs with genomic positions and provides back a tab delmited file of SNPs with genomic positions. The main function for this repository is the snp_sequence_puller_auto.sh, but it is dependent on snp_sequence_puller.sh. If snp_sequence_puller.sh and snp_sequence_puller_auto.sh are not found in the same directory, then snp_sequence_puller_auto.sh will throw an error message and fail to complete.   

### Requirments
There are several requirments of both functions. Below I will detail the required inputs of both functions.

#### snp_sequence_puller.sh
The snp_sequence_puller.sh function requires several inputs to properly function. Below is a discritpion of the flags.

1. **-f, --genome-file:**  
   - Input reference genome file. This should be a fasta (.fa) file and the fasta file should come with an index (.fai). If the index is not present, this code will create one. If you do not have writing privileges in the directory of your reference genome, then the function may fail due to not being able to write an index if one is not present.

2. **-p, --position:**  
   - Position in the reference genome. This should be an integer ranging from 0 to the end of the chromosome length.

3. **-l, --length:**  
   - Flanking sequence length. This is the length in base pairs around your SNP. Note, if you provide a length longer than the amount of base pairs on either side, you will be greeted with an error message.

4. **-c, --chromosome:**  
   - Chromosome in reference genome. This must match the chromosome names provided in the reference (e.g., Chr1A vs. 1A).

5. **-a, --alternate-allele:**  
   - Alternate allele of provided position. This must be a single nucleotide sequence. This must be an A, T, G, or C.

6. **-r, --reference-allele:**  
   - Reference allele of provided position. This must be a single nucleotide sequence. This must be an A, T, G, or C.

#### snp_sequence_puller_auto.sh
The snp_sequence_puller_auto.sh uses snp_sequence_puller.sh to recursivly perform the above task for every listed sequence in a tab delimited file. Below is a discription of the flags.

1. **-i, --input-file:**  
    - Input tab delmited file (e.g., input_file.txt)
    
2. **-o, --output-file:**  
    - Name of the tab delimited output file (e.g., output_file.txt)

3. **-l, --length:**
    - Length in bp (default is 200). This is the length of the flanking sequence you want around your SNP. 

4. **-r, --reference-geno:**  
    - Input reference genome file. This should be a fasta (.fa) file and the fasta file should come with an index (.fai). If the index is not present, this code will create one. If you do not have writing privileges in the directory of your reference genome, then the function may fail due to not being able to write an index if one is not present.

The input file (-i) must be a tab delmited file that has the following columns
1. chr (chromosome)
2. pos (position)
3. id (name of the SNP)
4. ref (reference allele)
5. alt (alternate allele)

To properly format, refer to the [example proivded](https://github.com/zjwinn/SNP-Sequence-Puller/blob/main/input_file_example.txt). Note, a sanity check will be performed to check if the reference provided and the "reference allele" provided agree. If this is not the case, the marker will be thrown out. 

### Usage
To use the shell scripts provided, call directly on them using the following code:
```bash
bash snp_sequence_puller_auto.sh
```
This will call on the script directly in the directory you are running in or it can be used to call on a directory outside the current working directory if specified. An example of how to run the function is provided below:
```bash
bash snp_sequence_puller_auto.sh \
   -i input_file_example.txt \
   -o example_output.txt \
   -l 200 \
   -r path_to_refrence_sequence.fa \
   -v > snp_sequence_puller_auto.sh.log
```
To run this example you will have to download the [RefSeqv2.1 Chinese Spring Wheat Genome Sequence](https://urgi.versailles.inrae.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v2.1/) and **change "path_to_refrence_sequence.fa" in the above script to the appropriate path**. If at any time you wish to display the help function, use the following:
```bash
bash snp_sequence_puller_auto.sh --help
```
### Package Requirments
This package requires the following to be installed:
1. bedtools - [link to page](https://bedtools.readthedocs.io/en/latest/index.html)
2. samtools/bcftools - [link to page](https://samtools.github.io/bcftools/)

Please make sure you have both programs installed properly prior to running to avoid errors.
