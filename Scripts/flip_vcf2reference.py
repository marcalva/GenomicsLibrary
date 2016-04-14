#!/usr/bin/python
#
# Force SNPs to be reference allele, when possible. Write out SNP IDs from forward, reverse, and unknown strand to separate file

import sys
import math

# First argument is vcf file
# Second argument is fasta file. Needs fasta index file with .fai extension
# Third argument is output vcf file
# Fourth argument is output SNP info file

class fai:
	def __init__(self, chr, length, offset, length_bases, line_bytes):
		self.chr = chr
		self.length = int(length)
		self.offset = int(offset)
		self.length_bases = int(length_bases)
		self.line_bytes = int(line_bytes)


def reverse(x):
	if x.upper() == "A":
		return "T"
	elif x.upper() == "C":
		return "G"
	elif x.upper() == "G":
		return "C"
	elif x.upper() == "T":
		return "A"
	else:
		return "NA"

vcf_fs = open(sys.argv[1], 'r')
in_fasta = open(sys.argv[2], 'r')
in_fasta_fai = open(sys.argv[2] + ".fai", 'r')

outvcf_fs = open(sys.argv[3], 'w')
outstat_fs = open(sys.argv[4], 'w')


fai_dictionary = {}
for line in in_fasta_fai:
	line_list = line.strip().split()
	thisfai = fai(line_list[0], line_list[1], line_list[2], line_list[3], line_list[4])
	fai_dictionary[line_list[0]] = thisfai


for line in vcf_fs:
	if line[0] == "#":
		outvcf_fs.write(line)
		continue
	
	line_list = line.strip().split()
	
	if len(line_list) < 10:
		print("vcf line has less than 10 elements")
		print(line)
		sys.exit(1)
	
	chr = line_list[0]
	bp = int(line_list[1])
	id = line_list[2]
	ref = line_list[3]
	alt = line_list[4]
	
	# Get reverse
	ref_reverse = reverse(ref)
	alt_reverse = reverse(alt)
	if (ref_reverse == "NA" or alt_reverse == "NA"):
		outstat_fs.write(id + "\tunknown\n")
		continue
	
	# Get position in fasta file
	thisfai = fai_dictionary[chr]
	offset = thisfai.offset
	length_bases = thisfai.length_bases
	length_bytes = thisfai.line_bytes
	bp0 = bp - 1
	line_number = math.floor(bp0 / length_bases)
	line_position = bp0 % length_bases
	snp_offset = offset + (line_number * length_bytes) + line_position
	
	# Get genome base
	in_fasta.seek(int(snp_offset))
	genome_base = in_fasta.read(1)
	genome_base_reverse = reverse(genome_base)
	
	if (genome_base == ref):
		outstat_fs.write(id + "\tforward\n")
		outvcf_fs.write(line)
	elif (genome_base == ref_reverse):
		outstat_fs.write(id + "\treverse\n")
		out_line = "\t".join([chr, str(bp), id, ref_reverse, alt_reverse] + line_list[5:])
		outvcf_fs.write(out_line + "\n")
#	elif (genome_base == alt):
#		outstat_fs.write(id + "\tforward/switch\n")
#	elif (genome_base == alt_reverse):
#		outstat_fs.write(id + "\treverse/switch\n")
#		out_line = "\t".join([chr, str(bp), id, ref_reverse, alt_reverse] + line_list[5:])
#		outvcf_fs.write(out_line + "\n")
	else:
		outstat_fs.write(id + "\tunknown\n")
	
vcf_fs.close()
outvcf_fs.close()
outstat_fs.close()
