#!/usr/bin/python
#
# Check allelic concordance between a given eQTL SNP and ASE SNPs.

import sys
import argparse
import re

#######################################
# Parse and set up argument parameters
#######################################


parser = argparse.ArgumentParser()
parser.add_argument("--vcf", "-v", help="Phased VCF file. If '-', then takes input from standard input.", required=True)
parser.add_argument("--eqtl", "-e", help="eQTL ID. Checks concordance against direction of this SNP.", required=True)
parser.add_argument("--ase", "-a", nargs='+', help="Set of ASE SNP IDs separated by spaces.", required=True)
parser.add_argument("--direction", "-d", help="Must be either '+' or '-'. Indicates direction of effect of alternate allele. + would indicate a higher expression level if the eQTL SNP has the alternate allele.", required=True)
parser.add_argument("--ase_flag", help="Flag in the sample field indicating the field number of the ASE counts. Defaults to 'AS'")
parser.add_argument("--min_cov", help="Minimum number of total reads at SNP to be used for calculations, defaults to 20"))

args = parser.parse_args()

if args.min_cov:
	minCov = args.min_cov
else:
	minCov = 20
	

if args.ase_flag:
	af = args.ase_flag
else:
	af = 'AS'

if args.direction == '+':
	ad = '+'
elif args.direction == '-':
	ad = '-'
else:
	print "Direction parameter must either be + or -"
	sys.exit(1)

eqtl = args.eqtl

ase = args.ase

if args.vcf == '-':
	vcfFs = sys.stdin
else:
	try:
		vcfFs = open(args.vcf, 'r')
	except:
		print "Could not open VCF file."



#######################################
# Store eQTL and ASE data
#######################################

# For eQTL we only need to store phased genotype data.
# For ASE data we need to store the phased count data (genotype is meaningless)

eqtlGeno = [] # List of tuples of genotype, coded 0 for ref and 1 for alt
hetIndex = []
asCounts = [] # list of list of tuples of allele counts

for line in vcfFs:
	if line[0] == "#":
		continue
	lineList = line.strip().split()
	
	if lineList[0] == eqtl:
		for i in range(9,len(lineList)):
			geno = lineList[i].strip().split(':')[0]
			hap = re.split('[/|]', geno)
			if hap[0] != hap[1]:
				continue
			hapTuple = (hap[0], hap[1])
			eqtlGeno.append(hapTuple)
			hetIndex.append(i)
	
	if lineList[0] in ase:
		thisCount = []
		genoField = lineList[8].strip().split(':')
		try:
			aCol = genoField.index(af)
		except:
			print "Could not find Allele Counts ID in Format Field."
		
		for i in hetIndex:
			counts = lineList[i].strip().split(':')[aCol]
			thisCount.append(counts[0], counts[1])
		
		asCounts.append(thisCount)

if len(eqtlGeno) != 1:
	print "Need one eQTL SNP for analysis"

if len(asCounts) < 1:
	print "No SNPs with allele counts matching eQTL"
	sys.exit(1)

if len(eqtlGeno) != len(asCounts[0]):
	print "Not the same number of samples at every SNP"

#######################################
# Calculate and output statistics
#######################################
# Calculate number/total and percent of ALL SNPs that match direction of eqtl
# For each SNP, calculate number/total and percent that match direction of eqtl


for fSnp in range(0,len(ase)):
	for sample in range(0, len(eqtlGeno)):
		hap1Counts
	