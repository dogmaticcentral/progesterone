#!/usr/bin/python3
# the single input file from endometrial microvascular endothelial cells from
# Job Dekker lab, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105710
# the (bed) file with TADs can be found here
# https://www.encodeproject.org/experiments/ENCSR551IPY/
# (under processed data)

import os
import subprocess
from statistics import mean, stdev
from math import floor, ceil

from linkto_python_modules.utils import *

#########################################
def main():

	gene_name = "Hand2"

	tadfile = "/storage/databases/tads/encode/ENCFF633ORE.bed"
	ucsc_gene_regions_dir = "/storage/databases/ucsc/gene_ranges/human/hg19"

	for prerequisite in [tadfile, ucsc_gene_regions_dir]:
		if os.path.exists(prerequisite): continue
		print(prerequisite, "not found")
		exit()

	chromosome, strand, gene_range = ucsc_gene_coords(gene_name, ucsc_gene_regions_dir)
	if not chromosome or not gene_range:
		exit()

	print ( "{} {} {}:{}-{}".format(gene_name, strand, chromosome, gene_range[0], gene_range[1]) )

	[start, end] = get_tad (tadfile, chromosome, gene_range)
	print ("TAD containing %s region: %s:%d-%d   length %d"%(gene_name, chromosome, start, end, end-start))


#########################################
if __name__ == '__main__':
	main()




