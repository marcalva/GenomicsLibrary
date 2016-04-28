#!/usr/bin/python
#
# Get GC content of genes based on BED intervals
# Input BED file is argument 1
# FASTA with index is argument 2
# Output file is argument 3

import sys
import gzip
import math

class fai:
	def __init__(self, chr, length, offset, length_bases, line_bytes):
		self.chr = chr
		self.length = int(length)
		self.offset = int(offset)
		self.length_bases = int(length_bases)
		self.line_bytes = int(line_bytes)


if len(sys.argv) != 3:
	print "Usage:"
	print "python get_gc_percentage.py BED_file FASTA_file OUT_file"
	print "BED file has coordinates for each feature (typically exons) in the 4th column"
	print "FASTA file must be indexed, and have the sequence for all features in the BED file"
	print "\tIndex by running 'samtools faidx <ref.fa>', outputting <ref.fa>.fai"
	print "Output file is tab delimited text file with feature name in first column and GC% in second column."
	sys.exit(1)

if sys.argv[1][-2:] == "gz":
	ifs = gzip.open(sys.argv[1], 'r')
else:
	ifs = open(sys.argv[1], 'r')

fasta = open(sys.argv[2], 'r')
fasta_fai = open(sys.argv[2] + ".fai", 'r')
fai_dictionary = {}
for line in fasta_fai:
	line_list = line.strip().split()
	thisfai = fai(line_list[0], line_list[1], line_list[2], line_list[3], line_list[4])
	fai_dictionary[line_list[0]] = thisfai
fasta_fai.close()

gene_gc_count = {}
genes_total_count = {}

for line in ifs:
	if line[0] == "#":
		continue
	
	line_list = line.strip().split()
	
	gene = line_list[3]
	chr = line_list[0]
	start = int(line_list[1])
	end = int(line_list[2])
	
	thisfai = fai_dictionary[chr]
	offset = thisfai.offset
	length_bases = thisfai.length_bases
	length_bytes = thisfai.line_bytes

	gc_count = 0
	total_count = 0
	# BED file is 0-based as is the FASTA file access.
	for bp in range(start, end):
		line_number = math.floor(bp / length_bases)
		line_position = bp % length_bases
		bp_offset = offset + (line_number * length_bytes) + line_position
		fasta.seek(int(bp_offset))
		nucleotide = fasta.read(1)
		# print "%s:%d %s" %(chr, bp, nucleotide)
		if nucleotide == "G" or nucleotide == "C":
			gc_count = gc_count + 1
		total_count = total_count + 1

	if gene not in gene_gc_count:
		gene_gc_count[gene] = gc_count
		genes_total_count[gene] = total_count
	else:
		gene_gc_count[gene] = gene_gc_count[gene] + gc_count
		genes_total_count[gene] = genes_total_count[gene] + total_count

ifs.close()
fasta.close()
ofs = open(sys.argv[3], 'w')
for g in gene_gc_count:
	gc_ratio = float(gene_gc_count[g]) / float(genes_total_count[g])
	outline = g + "\t" + str(gc_ratio) + "\n"
	ofs.write(outline)

ofs.close()
