#!/bin/bash

# Run test file 
bash ./Wheat-KASP-Designer.sh \
    --input-file 'NC13-20076xGA06493-13LE6_filt.vcf.gz' \
    --reference-geno  /90daydata/gbru_wheat2/zjwinn_project_directory/Ref/161010_Chinese_Spring_v1.0_pseudomolecules.fasta \
    --output-file 'test.txt' \
    --snps 'test_snps.txt' \
    --verbose \
    --debug
