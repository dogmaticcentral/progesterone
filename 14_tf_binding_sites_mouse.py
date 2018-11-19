#!/usr/bin/python3

# the analysis here is a bit smaller - only estrogen and progesterone receptor binding sites
# we try to produce similar table as in the case of human, using ChIPseq data for murine uterus


import sys, os

# pycharm recognizes this if it says .linkto
# however python3 does not like the dot
from linkto_python_modules.mysqldb import *
from linkto_python_modules.utils import *

def overlap(interval_list, qry_interval):

	ovlp = False
	for interval in interval_list:
		if qry_interval[1]<interval[0]: continue
		if interval[1]<qry_interval[0]: continue
		# some leeway coulf be left here one day ...
		ovlp = True
	return ovlp

def process_binding(filetype, data_dir, agonist_file, vehicle_file, chrom, region_start, region_end):

	chrom_column = 0

	if filetype=="pr":
		chrom = chrom.replace("chr","")
		chrom_column = 2

	elif filetype=="er":
		if type(chrom)==int:
			chrom = "chr%d"%chrom
		elif not 'chr' in chrom:
			chrom = 'chr'+chrom
		chrom_column = 0
	else:
		print("unrecognized file type:", filetype)
		exit()
		
	# progesterone
	agonist_binding_intervals = []
	inf = open("{}/{}".format(data_dir, agonist_file, "r"))
	for line in inf:
		if  filetype=="pr" and 'Sample' in line: continue
		if  filetype=="er" and 'track' in line: continue
		fields = line.rstrip().split("\t")
		if filetype=="pr" and len(fields)<5: continue
		if filetype=="er" and len(fields)<4: continue
		if fields[chrom_column] != str(chrom): continue
		[start,end] = [int(i) for i in fields[chrom_column+1:chrom_column+3]]
		if end<=region_start or region_end<=start: continue
		agonist_binding_intervals.append([start,end])
	inf.close()

	# vehicle
	vehicle_binding_intervals = []
	inf = open("{}/{}".format(data_dir, vehicle_file, "r"))
	for line in inf:
		if  filetype=="pr" and 'Sample' in line: continue
		if  filetype=="er" and 'track' in line: continue
		fields = line.rstrip().split("\t")
		if filetype=="pr" and len(fields)<5: continue
		if filetype=="er" and len(fields)<4: continue
		if fields[chrom_column] != str(chrom): continue
		[start,end] = [int(i) for i in fields[chrom_column+1:chrom_column+3]]
		if end<=region_start or region_end<=start: continue
		vehicle_binding_intervals.append([start,end])
	inf.close()

	for interval in agonist_binding_intervals:
		if overlap(vehicle_binding_intervals, interval):
			agonist_binding_intervals.remove(interval)
			continue

	return agonist_binding_intervals

#########################################
def main():

	assembly = "mm9"

	if len(sys.argv) < 2:
		print  ("usage: %s <gene_name> " % sys.argv[0])
		exit()

	gene_name = sys.argv[1]

	ucsc_gene_regions_dir = "/storage/databases/ucsc/gene_ranges/mouse/%s" % assembly

	pr_binding_dir = "/storage/databases/geo/GSE34927"
	pr_binding_P4_file = "GSM857546_2_PR_P4_s_1_aligned.tsv"
	pr_binding_vehicle = "GSM857545_1_PR_oil_s_4_aligned.tsv"

	er_binding_dir = "/storage/databases/geo/GSE36455"
	er_binding_E2_file = "GSM894054_WT_E2.ER_peaks.bed"
	er_binding_vehicle = "GSM894053_WT_vehicle.ER_peaks.bed"

	for dep in [ucsc_gene_regions_dir, pr_binding_dir, er_binding_dir,
				"{}/{}".format(pr_binding_dir,pr_binding_P4_file),
				"{}/{}".format(pr_binding_dir,pr_binding_vehicle),
				"{}/{}".format(er_binding_dir,er_binding_E2_file),
				"{}/{}".format(er_binding_dir,er_binding_vehicle)]:
		if not os.path.exists(dep):
			print(dep,"not found")
			exit()

	chrom, strand, gene_range = ucsc_gene_coords(gene_name, ucsc_gene_regions_dir)
	print ("%s position:"%gene_name, chrom, strand, gene_range)

	# human TAD that includes Hand2 is 1440kbp longrm -rf
	# so let's say we search in the region  1Mbp on each side od the gene region in mouse
	region_start = gene_range[0] -  1400000
	region_end   = gene_range[1] +  1400000
	pr_binding_regions = process_binding("pr", pr_binding_dir, pr_binding_P4_file, pr_binding_vehicle,
											chrom, region_start, region_end)
	er_binding_regions = process_binding("er", er_binding_dir, er_binding_E2_file, er_binding_vehicle,
											chrom, region_start, region_end)

	outf = open ("raw_data/%s_tfbs_%s.tsv"%(gene_name,assembly),"w")
	outf.write("\t".join(["% chrom", "chromStart", "chromEnd", "name",]) + "\n")
	for interval in pr_binding_regions:
		outf.write("\t".join([chrom, str(interval[0]), str(interval[1]), "PR"]) + "\n")
	for interval in er_binding_regions:
		outf.write("\t".join([chrom, str(interval[0]), str(interval[1]), "ESR1"]) + "\n")

	outf.close()

	return True


#########################################
########################################
if __name__ == '__main__':
	main()

