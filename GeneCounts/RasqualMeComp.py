#!/u/local/apps/python/3.4.3/bin/python3.4

from scipy import stats
import tabix

class snp:
	def __init__(self, gene, name, p_val, chr, bp):
		self.gene = gene
		self.name = name
		self.p_val = p_val
		self.chr = chr
		self.bp = bp

meFn = "metsim494.imputed.clean.cis.txt" # Put full path here
meFs = open(meFn, 'r')
meFs.readline() # Read header
meGenePvals = {}

for line in meFs:
	lineList = line.strip().split()
	gene = lineList[6]
	name = lineList[2]
	pval = float(lineList[10])
	chr = int(lineList[0])
	bp = int(lineList[1])
	thisSnp = snp(gene, name, p_val, chr, bp)
	if gene not in meGenPvals:
		meGenePvals[gene] = thisSnp

meFs.close()
oFn = ""
oFs = open(oFn, 'w')
oFs.write("GENE\tMePval\tRasqualPval\tRasqualAse\n")

geneListFn = "genes_run.txt"
geneListFs = open(geneListFn, 'r')
for line in geneListFs:
	lineList = line.strip().split()
	gene = (lineList[0])
 	rAse = 'NA'
	rPval = 'NA'
	rChr = ''
	rBp = ''
	rName = ''
	rFn = "results/" + gene + ".txt"
	rFs = open(rFn, 'r')
	for rLine in rFs:
		rLineList = rLine.strip().split()
		tstat = float(rLineList[10])
		pval = stats.chisqprob(tstat,1)
		if rPval == 'NA':
			rPval = pval
			rChr = int(rLineList[2])
			rBp = int(rLineList[3])
			rName = rLineList[1]
			if rLineList[18] == 0:
				rAse = '0'
			else:
				rAse = '1'
		else:
			if pval < rPval:
				rPval = pval
				rName = rLineList[1]
				rChr = int(rLineList[2])
				rBp = int(rLineList[3])
	rFs.close()
	
	thisSnp = meGenePvals[gene]
	
	if gene in meGenePvals:
		thisSnp = meGenePvals[gene]
		
		# get R2
		r2 = "NA"
		if thisSnp.name == rName:
			r2 = 1
		else:
			rGenoList = []
			meGenoList = []
			# Open VCF file
			
			vcfName = "../" + str(rChr) + "/chr" + str(rChr) + "_metsim494.vcf.gz"
			tVcf = tabix.open(vcfName)
			
			rVcf = tVcf.querys(rChr - 1, rBp, rBp + 1)
			for i in range(9, len(rVcf)):
				if rVcf[i].split(":")[0] == "0|0":
					rGenoList.append(0)
				elif rVcf[i].split(":")[0] == "1|1":
					rGenoList.append(2)
				else:
					rGenoList.append(1)
			
			meVcf = tVcf.querys(thisSnp.chr - 1, thisSnp.bp, thisSnp.bp + 1)
			for i in range(9, len(meVcf)):
				if meVcf[i].split(":")[0] == "0|0":
					meGenoList.append(0)
				elif meVcf[i].split(":")[0] == "1|1":
					meGenoList.append(2)
				else:
					meGenoList.append(1)
			
			r2 = stats.pearsonr(rVcf,rVcf)[0] ^ 2
			
		outLine = "\t".join([gene, thisSnp.name, str(thisSnp.chr), str(thisSnp.bp), str(thisSnp.p_val), rName, str(rChr), str(rBp), str(rPval), rAse, str(r2)]) + '\n'
		oFs.write(outLine)
	else:
		outLine = "\t".join([gene, "NA", "NA", "NA", ">0.1", rName, str(rChr), str(rBp), str(rPval), rAse, str(r2)]) + '\n'
		oFs.write(outLine)

geneListFs.close()
oFs.close()