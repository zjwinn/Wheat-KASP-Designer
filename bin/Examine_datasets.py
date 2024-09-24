# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 11:07:10 2019

@author: lurib
"""
###This scripts is to examine the struture of the wheat genome sequence
### downloaded from Esembleplants as 4Gb fasta file
################################################################## 
from Bio import SeqIO
total = 0
for seq in SeqIO.parse("C:\CropScience\python codes\SNP_Primer_Pipeline-master\Triticum_aestivum.IWGSC.dna.toplevel.fa","fasta"):
  total += len(seq)
print (total) #14547261565 nucleotides

###################################################################
# This load teh data en print the chr, first and last 40 nucleotides of each Chro 
from Bio import SeqIO
for record in SeqIO.parse("C:\CropScience\python codes\SNP_Primer_Pipeline-master\Triticum_aestivum.IWGSC.dna.toplevel.fa", "fasta"):
    start_seq = record.seq[:40] # first 10 letters
    end_seq = record.seq[-40:] # last 10 letters
    print(record.id + " " + start_seq + "..." + end_seq)

## This checks the version of Biopython
import Bio
print(Bio.__version__)

## This tells you info about your working directory
import os
cwd = os.getcwd()