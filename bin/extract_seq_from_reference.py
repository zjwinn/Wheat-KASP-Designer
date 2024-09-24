#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  get_flanking_for_variaiton_sites.py
#  
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#

from subprocess import call
import sys

# function to get the flanking sequences
def get_flanking(range_file, flanking_file):
	reference = "/Library/WebServer/Documents/blast/db/nucleotide/IWGSC_v2_ChrU.fa"
	cmd = "blastdbcmd -entry_batch " + range_file + " -db " + reference + " > " + flanking_file
	print ("Command to get the flanking sequences for each SNP\n", cmd)
	call(cmd, shell=True)

# read input file
def get_seq_name(infile):
	seq_name_list = []
	for line in open(infile):
		line = line.strip()
		if not line: # skip blank lines
			continue
		seq_name_list.append(line.replace("\t", " "))
	return seq_name_list


# function to extract sequences from a fasta file 
def get_fasta(infile, seq_name_list):
	fasta = {} # dictionary for alignment
	n = 0 # for sequence name list index
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if not line:
				continue
			if line.startswith(">"):
				sequence_name = line.split()[0].lstrip(">")
				if sequence_name in seq_name_list[n]: # in case seq name mismatch
					sequence_name = seq_name_list[n]
					n += 1
			else:
				fasta.setdefault(sequence_name, "")
				fasta[sequence_name] += line.rstrip()
	return fasta


def main(args):
	reference_list = ["/Library/WebServer/Documents/blast/db/nucleotide/IWGSC_v2_ChrU.fa", 
	"/Library/WebServer/Documents/blast/db/nucleotide/IWGSC_CSS_AB-TGAC_UCW_v1.fa",
	"/Library/WebServer/Documents/blast/db/nucleotide/161010_Chinese_Spring_v1.0_pseudomolecules.fasta",
	"/Volumes/DATA2/databases/ncbi/nucleotide/Triticum_aestivum.TGACv1.dna.toplevel.fa",
	"/Users/galaxy/blastdb/IWGSC_v1.1_HC_20170706_cds.fasta"]
	infile = args[1]
	outfile = args[2]
	reference = reference_list[int(args[3]) - 1] # reference 1 or 2
	# step 1: get the flanking sequences
	flanking_file = "temp_flanking_seq.fa"
	cmd = "blastdbcmd -entry_batch " + infile + " -db " + reference + " > " + flanking_file
	call(cmd, shell=True)
	# step 2: replace the fasta sequence names with more information
	seq_name_list = get_seq_name(infile)
	seq_fasta = get_fasta(flanking_file, seq_name_list)
	# step 3: write the organized fasta
	out = open(outfile, "w")
	for i in seq_name_list:
		out.write(">" + i + "\n")
		out.write(seq_fasta[i] + "\n")
	out.close()
	
	return 0

if __name__ == '__main__':
	import sys
	sys.exit(main(sys.argv))
