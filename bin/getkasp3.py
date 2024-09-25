#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  getkasp
#
#  Copyright 2016 Junli Zhang <zhjl86@gmail.com>
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

# getkasp3.py use files with names starting with "flanking*" from step 5 of the pipeline
# getkasp3.py now uses all the variation sites that can differ all homeologs

# change the input wildcard and the paths of primer3 and muscle accordingly.
# NOTES: the output primer pair includes the common primer and the primer with SNP A or T, you need to change the 3' nucleartide to get the primer for the other SNP.
# changes
# 2022-09-15: add a simple score for auto-selection and only use positions that can diff all homeologs (variations instead of variations2)

# Usage: getkasp3.py maximum_primer_Tm maximum_primer_size whether_to_pick_anyway (1 is yes, 0 is NO)
# example:  ../bin/getkasp3.py 63 25  0

### Imported
from subprocess import call
import getopt, sys, os, re
from glob import glob
import copy
#########################
## input
max_Tm = sys.argv[1] # max Tm, default 63, can be increased in case high GC region
max_size = sys.argv[2] # max primer size, default 25, can be increased in case low GC region
pick_anyway = sys.argv[3] # pick primer anyway even if it violates specific constrains

# get all the raw sequences
raw = glob("flanking_temp_marker*") # All file names start from "flanking"
raw.sort()

iupac = {"R": "AG", "Y": "TC", "S": "GC", "W": "AT", "K": "TG", "M": "AC"}

#from sys import platform
def get_software_path(base_path):
	if sys.platform.startswith('linux'): # linux
		primer3_path = base_path + "/primer3_core"
		muscle_path = base_path + "/muscle"
	elif sys.platform == "win32" or sys.platform == "cygwin": # Windows...
		primer3_path = base_path + "/primer3_core.exe"
		muscle_path = base_path + "/muscle.exe"
	elif sys.platform == "darwin": # MacOSX
		primer3_path = base_path + "/primer3_core_darwin64"
		muscle_path = base_path + "/muscle3.8.31_i86darwin64"
	return primer3_path, muscle_path

# function to get reverse complement
def ReverseComplement(seq):
	s1 = "BDHKMNRSVWYATGCbdhkmnrsvwyatgc"
	s2 = "VHDMKNYSBWRTACGvhdmknysbwrtacg"
	seq_dict = {s1[i]:s2[i] for i in range(len(s1))}
	return "".join([seq_dict[base] for base in reversed(seq)])

# classes
class Primers(object):
	"""A primer set designed by Primer3"""
	def __init__(self):
		self.name = ""
		self.start = 0
		self.end = 0
		self.length = 0
		self.tm = 0.0
		self.gc = 0.0
		self.anys = 0.0
		self.three = 0.0
		self.hairpin = 0.0
		self.end_stability = 0.0
		self.seq = ""
		self.difthreeall = "NO" # whether 3' site can differ all
		self.difnum = 0
		self.direction = ""
	def formatprimer(self):
		formatout = "\t".join(str(x) for x in [self.direction, self.start, self.end, self.difnum, self.difthreeall, self.length, self.tm, self.gc, self.anys, self.three, self.end_stability, self.hairpin, self.seq, ReverseComplement(self.seq)])
		return(formatout)

class PrimerPair(object):
	"""A pair of primers designed by Primer3"""
	def __init__(self):
		self.left = Primers()
		self.right = Primers()
		self.compl_any = "NA"
		self.compl_end = "NA"
		self.penalty = "NA"
		self.product_size = 0
		self.score = 0 # a simple score for auto-selection

# simple Tm calculator
def Tm(seq):
	t=0
	for a in seq:
		if a=='A' or a=='T':
			t=t+2
		if a=='C' or a=='G':
			t=t+4
	return t

# calculate GC content of a sequence
def Calc_GC(seq):
	t = 0.0 # float
	for a in seq:
		if a=='C' or a=='G':
			t += 1
	return t / len(seq) * 100

# function to find the segments with largest Tm
def FindLongestSubstring(s1, s2):
	longestStart = 0
	longestEnd = 0
	largestTm = 0
	start = 0
	gaps = [i for i, c in enumerate(s1) if c=='-' or s2[i]=='-']
	gaps.append(len(s1))
	for gap in gaps:
		end = gap
		tm = Tm(s1[start:end])
		if  tm > largestTm:
			longestStart = start
			longestEnd = end
			largestTm = tm
		start = gap + 1
	nL = len(s1[:longestStart].replace("-","")) # number of bases on the left of the longest segments
	nR = len(s1[longestEnd:].replace("-",""))  # number of bases on the right of the longest segments
	return [longestStart, longestEnd, nL, nR]

# alignment score
def score_pairwise(seq1, seq2, gapopen = -4.0, gapext = -1.0, match = 1.0, mismatch = -1.0):
	score = 0
	gap = False
	for i in range(len(seq1)):
		pair = (seq1[i], seq2[i])
		if not gap:
			if '-' in pair:
				gap = True
				score += gapopen
			elif seq1[i] == seq2[i]:
				score += match
			else:
				score += mismatch
		else:
			if '-' not in pair:
				gap = False
				if seq1[i] == seq2[i]:
					score += match
				else:
					score += mismatch
			else:
				score += gapext
	return score

# gap differences for an alignment
def gap_diff(seq1, seq2): # equal length
	ngap = 0
	for i in range(len(seq1)):
		pair = seq1[i] + seq2[i] # "AA", "A-", "-A" or "--"
		if pair.count("-") == 1:
			ngap += 1
	return ngap

# get the list of sequences in the homeolog groups for comparison with current primer
def get_homeo_seq(fasta, target, ids, align_left, align_right):
	s1 = fasta[target] # primer sequence in the template with gaps
	seq2comp = [] # final sequence to compare for each other homeolog
	for k in ids:
		s2 = fasta[k]
		targetSeq = s1[align_left:(align_right + 1)]
		homeoSeq = s2[align_left:(align_right + 1)]
		score1 = score_pairwise(targetSeq, homeoSeq) # score in multiple alignment
		#print "Targetseq ", targetSeq
		#print "homeoSeq  ", homeoSeq
		# Get the sequences for comparison
		indexL, indexR, nL, nR = FindLongestSubstring(targetSeq, homeoSeq)
		indexL += align_left
		indexR += align_left
		seqL = s2[:indexL].replace("-","")
		seqR = s2[indexR:].replace("-","")
		#print "indexL, indexR, nL, nR ", indexL, indexR, nL, nR
		#print "s2[indexL:indexR] ", s2[indexL:indexR]
		if len(seqL) < nL: # in case it does not have enough bases
			seqL = "-" * (nL - len(seqL)) + seqL
		if len(seqR) < nR:
			seqR = seqR + "-" * (nR - len(seqR))
		seqk = seqL[::-1][:nL][::-1] + s2[indexL:indexR] + seqR[:nR]
		#print "primer   :", targetSeq.replace("-","")
		#print "seqk     :", seqk
		score2 = score_pairwise(targetSeq.replace("-",""), seqk)
		# if there are more than 3 gaps, the Tm usually will be 10 C lower than the perfect match
		# so just use gap shift 
		if score1 > score2 and gap_diff(targetSeq, homeoSeq) < 4:
			# print "homeoSeq but remove all the gaps"
			# print "targetSeq:", targetSeq
			# print "homeoSeq :", homeoSeq
			# print "seqk     :", seqk
			# print "primer   :", targetSeq.replace("-","")
			seqk = "".join([homeoSeq[i] for i, c in enumerate(targetSeq) if c!='-'])
		seq2comp.append(seqk)
		#print k, "\t", seqk
	return seq2comp

# function to count mismtaches
def mismatchn (s1, s2):
	return sum(c1!=c2 for c1,c2 in zip(s1,s2))

# function to extract sequences from a fasta file 
def get_fasta(infile):
	fasta = {} # dictionary for alignment
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if line: # skip blank lines
				if line.startswith(">"):
					sequence_name = line.lstrip("> ").split()[0] # left strip > or space, so " > abc edf" will be "abc edf", then split by space to get "abc"
					fasta[sequence_name] = ""
				else:
					fasta[sequence_name] += line.replace(" ", "") # remove spaces in case
	return fasta

# in case multiple hit in the same chromosome in the psudomolecule
# the input fasta file are in order from the blast file
def get_fasta2(infile, target_chrom):
	fasta = {} # dictionary for alignment
	target = ""
	non_target_list = []
	n = 0 # add a number in the chromosome name
	with open(infile) as file_one:
		for line in file_one:
			line = line.strip()
			if line.startswith(">"):
				sequence_name = line.split()[0].lstrip(">")
				sequence_name += "-" + str(n)
				n += 1
				if not target and target_chrom in sequence_name:
					target = sequence_name
				else:
					non_target_list.append(sequence_name)
			else:
				fasta.setdefault(sequence_name, "")
				fasta[sequence_name] += line.rstrip()
			
	return fasta, target, non_target_list

# parse primer3 output
def parse_primer3output(primer3output, primerpair_to_return):
	regex = "012345"[:primerpair_to_return]
	regex = "([" + regex + "])"
	primerpairs= {}
	with open(primer3output) as infile:
		for line in infile:
			line = line.strip()
			if "SEQUENCE_ID" in line:
				seqid = line.split("=")[1]
				for i in range(0, primerpair_to_return):
					primerpairs[seqid + "-" + str(i)] = PrimerPair()
				continue
			elif re.search("PRIMER_PAIR_" + regex + "_PENALTY", line):
				mm = re.search("PRIMER_PAIR_" + regex + "_PENALTY", line)
				primerpairs[seqid + "-" + mm.group(1)].penalty = line.split("=")[1]
				continue
			elif re.search( "PRIMER_PAIR_" + regex + "_COMPL_ANY", line):
				mm = re.search( "PRIMER_PAIR_" + regex + "_COMPL_ANY", line)
				primerpairs[seqid + "-" + mm.group(1)].compl_any = line.split("=")[1]
				continue
			elif re.search( "PRIMER_PAIR_" + regex + "_COMPL_END", line):
				mm = re.search( "PRIMER_PAIR_" + regex + "_COMPL_END", line)
				primerpairs[seqid + "-" + mm.group(1)].compl_end = line.split("=")[1]
			elif re.search( "PRIMER_PAIR_" + regex + "_PRODUCT_SIZE", line):
				mm = re.search( "PRIMER_PAIR_" + regex + "_PRODUCT_SIZE", line)
				primerpairs[seqid + "-" + mm.group(1)].product_size = int(line.split("=")[1])
				continue
			elif re.search( "PRIMER_LEFT_" + regex + "_SEQUENCE", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "_SEQUENCE", line)
				primerpairs[seqid + "-" + mm.group(1)].left.seq = line.split("=")[1]
				continue
			elif re.search( "PRIMER_LEFT_" + regex + "=", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "=", line)
				primerpairs[seqid + "-" + mm.group(1)].left.start = int(line.split("=")[1].split(",")[0])
				primerpairs[seqid + "-" + mm.group(1)].left.length = int(line.split("=")[1].split(",")[1])
				primerpairs[seqid + "-" + mm.group(1)].left.end = primerpairs[seqid + "-" + mm.group(1)].left.start + primerpairs[seqid + "-" + mm.group(1)].left.length - 1
				continue
			elif re.search( "PRIMER_LEFT_" + regex + "_TM", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "_TM", line)
				primerpairs[seqid + "-" + mm.group(1)].left.tm = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_LEFT_" + regex + "_GC_PERCENT", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "_GC_PERCENT", line)
				primerpairs[seqid + "-" + mm.group(1)].left.gc = float(line.split("=")[1])
			elif re.search( "PRIMER_LEFT_" + regex + "_SELF_ANY_TH", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "_SELF_ANY_TH", line)
				primerpairs[seqid + "-" + mm.group(1)].left.anys = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_LEFT_" + regex + "_SELF_END_TH", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "_SELF_END_TH", line)
				primerpairs[seqid + "-" + mm.group(1)].left.three = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_LEFT_" + regex + "_HAIRPIN_TH", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "_HAIRPIN_TH", line)
				primerpairs[seqid + "-" + mm.group(1)].left.hairpin = float(line.split("=")[1])
			elif re.search( "PRIMER_LEFT_" + regex + "_END_STABILITY", line):
				mm = re.search( "PRIMER_LEFT_" + regex + "_END_STABILITY", line)
				primerpairs[seqid + "-" + mm.group(1)].left.end_stability = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_RIGHT_" + regex + "_SEQUENCE", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "_SEQUENCE", line)
				primerpairs[seqid + "-" + mm.group(1)].right.seq = line.split("=")[1]
			elif re.search( "PRIMER_RIGHT_" + regex + "=", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "=", line)
				primerpairs[seqid + "-" + mm.group(1)].right.start = int(line.split("=")[1].split(",")[0])
				primerpairs[seqid + "-" + mm.group(1)].right.length = int(line.split("=")[1].split(",")[1])
				primerpairs[seqid + "-" + mm.group(1)].right.end = primerpairs[seqid + "-" + mm.group(1)].right.start - primerpairs[seqid + "-" + mm.group(1)].right.length + 1
				continue
			elif re.search( "PRIMER_RIGHT_" + regex + "_TM", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "_TM", line)
				primerpairs[seqid + "-" + mm.group(1)].right.tm = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_RIGHT_" + regex + "_GC_PERCENT", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "_GC_PERCENT", line)
				primerpairs[seqid + "-" + mm.group(1)].right.gc = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_RIGHT_" + regex + "_SELF_ANY_TH", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "_SELF_ANY_TH", line)
				primerpairs[seqid + "-" + mm.group(1)].right.anys = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_RIGHT_" + regex + "_SELF_END_TH", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "_SELF_END_TH", line)
				primerpairs[seqid + "-" + mm.group(1)].right.three = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_RIGHT_" + regex + "_HAIRPIN_TH", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "_HAIRPIN_TH", line)
				primerpairs[seqid + "-" + mm.group(1)].right.hairpin = float(line.split("=")[1])
				continue
			elif re.search( "PRIMER_RIGHT_" + regex + "_END_STABILITY", line):
				mm = re.search( "PRIMER_RIGHT_" + regex + "_END_STABILITY", line)
				primerpairs[seqid + "-" + mm.group(1)].right.end_stability = float(line.split("=")[1])
				continue
	return primerpairs

# function to find primer sequence variation site and highligh them in primer sequences
def format_primer_seq(primer, variation): # input is a primer object and variation list
	if primer.start < primer.end:
		start = primer.start
		end = primer.end
		seq = primer.seq.lower()
		#primer_range = range(primer.start - 1, primer.end)
	else:
		start = primer.end
		end = primer.start
		seq = ReverseComplement(primer.seq.lower())
		#primer_range = range(primer.end - 1, primer.start)
	
	primer_range = list(range(start - 1, end))
	var_sites = set(variation).intersection(primer_range)
	var_sites_relative = [i - start + 1 for i in var_sites]
	#seq = primer.seq.lower()
	#seq = primer.seq
	for i in var_sites_relative:
		seq = seq[:i] + seq[i].upper() + seq[i+1:]
	if primer.start < primer.end:
		primer.seq = seq
	else:
		primer.seq = ReverseComplement(seq)
	primer.difnum = len(var_sites)
	return primer


def kasp(seqfile):
	#flanking_temp_marker_IWB1855_7A_R_251.fa
	#snpname, chrom, allele, pos =re.split("_|\.", seqfile)[3:7]
	#snpname, chrom, allele, pos =re.split("_", seqfile[:-7])[3:7] # [:-7] remove .txt.fa in the file
	info = re.split("_", seqfile[:-7])
	snpname = info[3]
	chrom = "_".join(info[4:-2])
	allele = info[-2]
	pos = info[-1]
	#print "Pos ", pos
	snp_site = int(pos) - 1 # 0-based
	getkasp_path = os.path.dirname(os.path.realpath(__file__))
	global_setting_file = getkasp_path + "/global_settings.txt"
	directory = "KASP_output"
	if not os.path.exists(directory):
		os.makedirs(directory)
	out = directory + "/selected_KASP_primers_" + snpname + ".txt"
	product_min = 50
	product_max = 250
	alt_allele = iupac[allele][0] # choose A or T
	SNP_A, SNP_B = iupac[allele] # SNP 2 alleles
	#print "SNP_A, SNP_B ", SNP_A, SNP_B
	# software path
	primer3_path, muscle_path = get_software_path(getkasp_path)
	
	# get target and ids and rename fasta seq names
	fasta_raw, target, ids = get_fasta2(seqfile, chrom) # target is the target chromosome, ids are other non-target chromosome name list
	print(("target ", target))
	print(("others ", ids))
	# write the renamed fasta seq to a file
	seqfile2 = "renamed_" + seqfile
	out_temp = open(seqfile2, "w")
	for k, v in fasta_raw.items():
		out_temp.write(">" + k + "\n" + v + "\n")
	out_temp.close()
	
	seq_template = fasta_raw[target]
	variation = [] # variation sites that can differ ALL
	variation2 = [] # variation sites that can differ at least 2 homeologs
	
	# STEP 0: create alignment file and primer3output file
	RawAlignFile = "alignment_raw_" + snpname + ".fa"
	alignmentcmd = muscle_path + " -in " + seqfile2 + " -out " + RawAlignFile + " -quiet"
	print(("Alignment command: ", alignmentcmd))
	call(alignmentcmd, shell=True)
	settings_common = "PRIMER_TASK=generic" + "\n" + \
		"SEQUENCE_TEMPLATE=" + seq_template + "\n" + \
		"PRIMER_PRODUCT_SIZE_RANGE=50-100 100-150 150-250" + "\n" + \
		"PRIMER_THERMODYNAMIC_PARAMETERS_PATH=" + getkasp_path + "/primer3_config/"  + "\n" + \
		"PRIMER_MAX_SIZE=" + max_size + "\n" + \
		"PRIMER_MIN_TM=57.0" + "\n" + \
		"PRIMER_OPT_TM=60.0" + "\n" + \
		"PRIMER_MAX_TM=" + max_Tm + "\n" + \
		"PRIMER_PAIR_MAX_DIFF_TM=6.0" + "\n" + \
		"PRIMER_FIRST_BASE_INDEX=1" + "\n" + \
		"PRIMER_LIBERAL_BASE=1" + "\n" + \
		"PRIMER_NUM_RETURN=5"  + "\n" + \
		"PRIMER_EXPLAIN_FLAG=1"  + "\n" + \
		"PRIMER_PICK_ANYWAY=" + pick_anyway + "\n"

	###############################
	if not len(ids): # if there is no homeologs found, such as when ploidy is 1
		# loop to write primer3 input for each variation site
		# primer3 inputfile
		primer3input = directory + "/primer3.input." + snpname
		p3input = open(primer3input, 'w')
		# because A and T give lower Tm, so use them as template
		if alt_allele in "ATat":
			seq_template = seq_template[:snp_site] +  alt_allele + seq_template[snp_site + 1:]
		
		settings = settings_common + \
		"SEQUENCE_ID=" + snpname + "-left\n" + \
		"SEQUENCE_FORCE_LEFT_END=" + str(snp_site + 1) + "\n" + \
		"="	
		p3input.write(settings + "\n")
		
		settings = settings_common + \
		"SEQUENCE_ID=" + snpname + "-right\n" + \
		"SEQUENCE_FORCE_RIGHT_END=" + str(snp_site + 1) + "\n" + \
		"="	
		p3input.write(settings + "\n")
		p3input.close()
		# primer3 output file
		primer3output = directory + "/primer3.output." + snpname
		p3cmd = primer3_path + " -default_version=2 -output=" + primer3output + " -p3_settings_file=" + global_setting_file + " " + primer3input
		print(("Primer3 command 1st time: ", p3cmd))
		call(p3cmd, shell=True)
		primerpairs = parse_primer3output(primer3output, 5)
		# Get primer list for blast
		#primer_for_blast = {}
		final_primers = {}
		nL = 0 # left primer count
		nR = 0 # right primer count
		for i, pp in primerpairs.items():
			if pp.product_size != 0:
				pl = pp.left
				pr = pp.right
				pp.left = pl
				pp.right = pr
				final_primers[i] = pp
	else: # if there are homeolog sequences	
		########################
		# read alignment file
		fasta = get_fasta(RawAlignFile)
		# get the variaiton site among sequences
		print(("The target: ", target))
		print(("The other groups: ", ids))

		alignlen = len(fasta[target])
		print(("Alignment length: ", alignlen))
		
		# get the target ID template base coordinate in the alignment
		t2a = {} # template to alignment
		a2t = {}
		ngap = 0 # gaps
		for i in range(alignlen):
			if fasta[target][i] == "-":
				ngap += 1
				continue
			t2a[i - ngap] = i
			a2t[i] = i - ngap

		print(("last key of t2a", i - ngap))
		print(("last key of a2t", i))
		
		seq_template = fasta[target].replace("-","") # remove all gaps

		variation = [] # variation sites that can differ ALL
		variation2 = [] # variation sites that can differ at least 2 homeologs
		
		# calculate gap number on the left and right
		gap_left_target = len(fasta[target]) - len(fasta[target].lstrip('-'))
		gap_left = max([len(v) - len(v.lstrip('-')) for k, v in fasta.items()])
		gap_right = min([len(v.rstrip('-')) for k, v in fasta.items()])
		print(("gap_left_target, gap_left and gap_right: ", gap_left_target, gap_left, gap_right))
		
		diffarray = {} # a list of 0 or 1: the same as or different from the site in each sequences of ids
		#for i in range(alignlen):
		for i in range(gap_left, gap_right): # exclude 20 bases on each side
			b1 = fasta[target][i]
			if b1 == "-":  # target non-gap base
				ngap += 1
				continue
			pos_template = a2t[i] # position in the target template (no gaps)
			if pos_template < 20 or pos_template > len(seq_template) - 20:
				continue
			nd = 0 # number of difference
			da = [0] * len(ids) # differ array
			m = 0 # counter of next loop
			if pos_template < snp_site:
				align_left = t2a[pos_template - 19] # 20 bp left of current position
				align_right = i
			else:
				align_left = i # 20 bp left of current position
				align_right = t2a[pos_template + 19]
			seq2comp = get_homeo_seq(fasta, target, ids, align_left, align_right) # list of sequences for comparison
			for k in seq2comp:
				if pos_template < snp_site:
					b2 = k[-1] # homeolog non-gap bas
				else:
					b2 = k[0]
				if b1 != b2:
					nd += 1
					da[m] = 1 # m sequence has variation from target
				m += 1

			# for each site pos_template
			diffarray[pos_template] = da
			if nd == len(ids): # different from all other sequences
				if pos_template not in variation: # in case multiple gaps
					variation.append(pos_template)
			if nd > 0: # different from at least 1 other sequence
				if pos_template not in variation2: # in case multiple gaps
					variation2.append(pos_template)
		
		print(("Sites that can differ all\n", variation))
		print(("\nSites that can differ at least 1\n", variation2))
		#print "\nKeys of diffarray: ", diffarray.keys()
		#############
		# loop to write primer3 input for each variation site
		# primer3 inputfile
		primer3input = directory + "/primer3.input." + snpname
		p3input = open(primer3input, 'w')
		# because A and T give lower Tm, so use them as template
		if alt_allele in "ATat":
			seq_template = seq_template[:snp_site] +  alt_allele + seq_template[snp_site + 1:]

		for i in variation: # use variation2 if you want semi-specific primers.
			if i == snp_site:
				continue
			elif i < snp_site:
				left_end = i
				right_end = snp_site
			else:
				left_end = snp_site
				right_end = i
			if right_end - left_end > product_max - 35: # suppose both primers are 18 bp
				continue
			settings = settings_common + \
			"SEQUENCE_ID=" + snpname + "-" + str(i+1) + "\n" + \
			"SEQUENCE_FORCE_LEFT_END=" + str(left_end + 1) + "\n" + \
			"SEQUENCE_FORCE_RIGHT_END=" + str(right_end + 1) + "\n" + \
			"="
			p3input.write(settings + "\n")

		p3input.close()

		# primer3 output file
		primer3output = directory + "/primer3.output." + snpname
		p3cmd = primer3_path + " -default_version=2 -output=" + primer3output + " -p3_settings_file=" + global_setting_file + " " + primer3input
		print(("Primer3 command 1st time: ", p3cmd))
		call(p3cmd, shell=True)
		primerpairs = parse_primer3output(primer3output, 1)
		print(("primerpairs length ", len(primerpairs)))
		####################################
		# Get primer list for blast
		#primer_for_blast = {}
		final_primers = {} # final primers for output
		nL = 0 # left primer count
		nR = 0 # right primer count
		for i, pp in primerpairs.items():
			dif3all = 0 # whether common primer can diff all in the 3' end
			varsite = int(i.split("-")[-2]) - 1 # variation site
			#print "varsite", varsite
			if pp.product_size != 0:
				pl = pp.left
				pr = pp.right
				# check whether 3' can differ all
				if varsite in variation:
					pl.difthreeall = "YES"
					pr.difthreeall = "YES"
					dif3all = 1
				if varsite < snp_site:
					pc = pl # pc is the common primer
					# print "pc = pl"
					# rr: range to check; only check 10 bases from 3' end
					rr = list(range(max(pc.end - 10,gap_left), pc.end)) # pc.end is 1 based, so change to 0 based.
				else:
					pc = pr
					# print "pc = pr"
					rr = list(range(pc.end -1, min(pc.end + 9, len(seq_template) - 20))) # rr should be within the keys of diffarray, which is from gap_left to gap_right
				# print "gap_left ", gap_left
				# print "len(seq_template) ", len(seq_template)
				# print "pc.end ", pc.end
				# print "rr ", rr

				# calculate a simple score for selection: only consider product size, 3'diff all, and Tm diff, diff number
				pp.score = dif3all*5.0 + 150.0/pp.product_size + pc.difnum/10.0 - abs(pl.tm - pr.tm)/10.0

				# sum of all the variation in each site
				aa = [sum(x) for x in zip(*(diffarray[k] for k in rr))]
				#print "aa ", aa
				if min(aa) > 0: # if common primer can differ all
					pp.left = pl
					pp.right = pr
					final_primers[i] = pp
	#################################################	
	# write to file
	outfile = open(out, 'w')
	outfile.write("index\tproduct_size\ttype\tstart\tend\tvariation number\t3'diffall\tlength\tTm\tGCcontent\tany\t3'\tend_stability\thairpin\tprimer_seq\tReverseComplement\tpenalty\tcompl_any\tcompl_end\tscore\n")			
	# write output file
	for i, pp in final_primers.items():
		pl = format_primer_seq(pp.left, variation)
		pr = format_primer_seq(pp.right, variation)
		pl.direction = "LEFT"
		pr.direction = "RIGHT"
		if pl.end == snp_site + 1:
			#print "pl.seq, SNP_A, SNP_B, ReverseComplement(SNP_A), ReverseComplement(SNP_B)"
			#print pl.seq, SNP_A, SNP_B, ReverseComplement(SNP_A), ReverseComplement(SNP_B)
			pA = copy.deepcopy(pl)
			pB = copy.deepcopy(pl)
			pA.seq = pA.seq[:-1] + SNP_A
			pB.seq = pB.seq[:-1] + SNP_B
			pC = pr
		else:
			#print "pr.seq, SNP_A, SNP_B, ReverseComplement(SNP_A), ReverseComplement(SNP_B)"
			#print pr.seq, SNP_A, SNP_B, ReverseComplement(SNP_A), ReverseComplement(SNP_B)
			pA = copy.deepcopy(pr)
			pB = copy.deepcopy(pr)
			pA.seq = pA.seq[:-1] + ReverseComplement(SNP_A)
			pB.seq = pB.seq[:-1] + ReverseComplement(SNP_B)
			pC = pl
		outfile.write("\t".join([i + "-Allele-" + SNP_A, str(pp.product_size), pA.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, str(pp.score)]) + "\n")
		outfile.write("\t".join([i + "-Allele-" + SNP_B, str(pp.product_size), pB.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, str(pp.score)]) + "\n")
		outfile.write("\t".join([i + "-Common", str(pp.product_size),   pC.formatprimer(), pp.penalty, pp.compl_any, pp.compl_end, str(pp.score)]) + "\n")

	outfile.write("\n\nSites that can differ all for " + snpname + "\n")
	outfile.write(", ".join([str(x + 1) for x in variation])) # change to 1 based
	outfile.write("\n\n\n")
	outfile.close()

# loop for all snp sequence files

for ff in raw:
	kasp(ff)

