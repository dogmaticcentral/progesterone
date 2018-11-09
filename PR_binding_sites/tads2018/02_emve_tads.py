#!/usr/bin/python3
# the single input file from endometrial microvascular endothelial cells from
# JOb Dekker lab, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105710
# the (bed) file with TADs can be found here
# https://www.encodeproject.org/experiments/ENCSR551IPY/
# (under processed data)

import os
import subprocess
from statistics import mean, stdev
from math import floor, ceil

#########################################
def  process_ucsc_ret(gene_name, ucsc_gene_regions_dir, ret):
	lines = []
	for line in ret.split("\n"):
		fields = line.split("\t")
		[infile, refseq_id] = fields[0].split(":")
		lines.append(fields[1:])
	if len(lines)==0:
		print("no entry for %s found in %s " % (gene_name, ucsc_gene_regions_dir))
		return None, None
	if len(lines)==2:
		print("more than one entry found for %s found in %s " % (gene_name, ucsc_gene_regions_dir))
		return None, None
	# we assume certain format in the file name, containting the chromosome number: e.g. chr18.csv
	chromosome = infile.split("/")[-1].replace(".csv","")
	[name, strand, txStart, txEnd] = lines[0]
	return chromosome, [int(txStart), int(txEnd)]

#########################################
def main():

	gene_name = "Hand2"

	tadfile = "/storage/databases/tads/encode/ENCFF633ORE.bed"
	ucsc_gene_regions_dir = "/storage/databases/ucsc/gene_ranges/human/hg19"

	for prerequisite in [tadfile, ucsc_gene_regions_dir]:
		if os.path.exists(prerequisite): continue
		print(prerequisite, "not found")
		exit()

	cmd = "grep -i %s %s/*" % (gene_name, ucsc_gene_regions_dir)
	ret =  subprocess.check_output(cmd, shell=True).decode('utf-8').rstrip()
	if not ret or len(ret)==0:
		print ("no entry for %s found in %s " % (gene_name, ucsc_gene_regions_dir))
		exit()

	chromosome, gene_range = process_ucsc_ret(gene_name, ucsc_gene_regions_dir,ret)
	if not chromosome or not gene_range:
		exit()

	print (gene_name, chromosome, gene_range)

	tads = {}
	inf = open(tadfile,"r")
	for line in inf:
		[chr, start, end] = line.rstrip().split()[:3]
		if not chr in tads: tads[chr] = []
		tads[chr].append([int(start), int(end)])


	for start,end in tads[chromosome]:
		if start<=gene_range[0]<=end or start<=gene_range[1]<=end:
			print (start, end, end-start)

#########################################
if __name__ == '__main__':
	main()




