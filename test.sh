#!/bin/bash

# Run test file 
bash ./Wheat-KASP-Designer.sh \
    --input-file 'NC13-20076xGA06493-13LE6_filt.vcf.gz' \
    --reference-geno '/mnt/c/Users/zachary.winn/Downloads/iwgsc_refseqv2.1_assembly.fa' \
    --output-file 'test.txt' \
    --snps 'test_snps.txt' \
    --verbose \
    --debug
