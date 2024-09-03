# Wheat KASP Designer

## Introduction
This repository contains a submodules of the popular [primer3](https://github.com/primer3-org/primer3) Git repository. In this repository, I (@zjwinn) will be writing a pipeline which may be used to take a properly formatted variant calling format (VCF) file and a reference genome (e.g., [Chinese Spring RefSeq V2.1](https://onlinelibrary.wiley.com/doi/10.1111/tpj.15289))output potential Kompetative allele specific primer ([KASP](https://www.biosearchtech.com/support/education/kasp-genotyping-reagents/how-does-kasp-work)) designs for biallelic single nucleotide polymoprhisms (SNPs). The purpose of this repository is to create an exicutable batch file (.sh) which can be called upon to design KASP automatically for mapping population or genome-wide marker panels. 

## Primer3 Information

Design PCR primers from DNA sequence. Widely used (190k Google hits for "primer3").
From mispriming libraries to sequence quality data to the generation of internal
oligos, primer3 does it. C&perl.

### Installing

```{bash}
`sudo apt-get install -y build-essential g++ cmake git-all`

`git clone https://github.com/primer3-org/primer3.git primer3`

`cd primer3/src`

`make`

`make test`
```

### Run Primer3

```{bash}
`./primer3_core ../example`
````

### Read the complete Primer3 manual

[Primer3 Manual](http://primer3.org/manual.html)
