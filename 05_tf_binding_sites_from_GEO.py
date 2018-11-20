#!/usr/bin/python3

# the main problem with the local (i.e. downloaded from GEO)
# is the arbitrariness of the format

import sys, os

# pycharm recognizes this if it says .linkto
# however python3 does not like the dot
from utils.mysqldb import *
from utils.utils import *
from utils.CrossMap import *

def overlap(interval_list, qry_interval):
	ovlp = False
	for interval in interval_list:
		if qry_interval[1]<interval[0]: continue
		if interval[1]<qry_interval[0]: continue
		# some leeway could be left here one day ...
		ovlp = True
	return ovlp


def read_bed(infile, region_chrom, region_start, region_end):
	intervals = []
	inf = open(infile, "r")
	for line in inf:
		fields = line.rstrip().split("\t")
		if len(fields)<3: continue
		chrom = fields[0].replace('chr','')
		if chrom != region_chrom: continue
		try:
			[start,end] = [int(i) for i in fields[1:3]]
		except:
			continue
		if (region_start and region_end) and (end<=region_start or region_end<=start): continue
		intervals.append([start,end])
	inf.close()
	return intervals

#############
def read_binding_intervals(data_dir, agonist_file, vehicle_file, chrom, region_start, region_end):
	# agonist
	infile = "{}/{}".format(data_dir, agonist_file)
	agonist_binding_intervals = read_bed(infile, chrom, region_start, region_end)

	if not vehicle_file: return agonist_binding_intervals

	# if we have control file, subtract regions that popped up
	# with vehicle only:
	infile = "{}/{}".format(data_dir, vehicle_file)
	vehicle_binding_intervals = read_bed(infile, chrom, region_start, region_end)

	for interval in agonist_binding_intervals:
		if overlap(vehicle_binding_intervals, interval):
			agonist_binding_intervals.remove(interval)
			continue

	return agonist_binding_intervals

#########################################
def process_tf_binding(outdir, geodir, tadfile, input_line):

	[species, gene_name, tf_name, geo_id, agonist_file, ctrl_file, input_assembly] = input_line.rstrip().split("\t")

	#############################
	# input checking
	chain_file=None
	ref_assembly = None
	if species=="human":
		ref_assembly="hg19"
	elif species=="mouse":
		ref_assembly="mm9"
	else:
		print("organism \"%s\" not recognized"%species)
		exit()

	if input_assembly != ref_assembly:  # we'll need  tools to translate
		chain_file="/storage/databases/liftover/{}To{}.over.chain".format(input_assembly, ref_assembly.capitalize())


	ucsc_gene_regions_dir = "/storage/databases/ucsc/gene_ranges/%s/%s" % (species, ref_assembly)

	if ctrl_file.replace(" ","") == "": ctrl_file = None
	data_dir = "{}/{}".format(geodir, geo_id)
	dependencies = [ucsc_gene_regions_dir, data_dir, "{}/{}".format(data_dir,agonist_file)]
	if ctrl_file: dependencies += ["{}/{}".format(data_dir,ctrl_file)]
	if chain_file: dependencies += [chain_file]
	for dep in dependencies:
		if not os.path.exists(dep):
			print(dep,"not found")
			exit()

	#############################
	# fetching additional info
	if 'chr' in gene_name:
		chrom = gene_name
		gene_name = None
		region_start, region_end = None, None
	else:
		chrom, strand, gene_range = ucsc_gene_coords(gene_name, ucsc_gene_regions_dir)
		print ("%s position:"%gene_name, chrom, strand, gene_range)
		if species=="human":
			[region_start, region_end] = get_tad (tadfile, chrom, gene_range)
		else:
			# human TAD that includes Hand2 is 1440kbp longrm -rf
			# so let's say we search in the region  1Mbp on each side od the gene region in mouse
			region_start = gene_range[0] - 1400000
			region_end   = gene_range[1] + 1400000

	#############################
	# business
	binding_regions = read_binding_intervals(data_dir, agonist_file, ctrl_file, chrom.replace('chr', ''), region_start, region_end)

	#############################
	# output
	outfile = "%s/%s_%s_tfbs_%s.tsv"%(outdir, gene_name if gene_name else chrom, tf_name, input_assembly)
	outf = open (outfile,"w")
	outf.write("\t".join(["% chrom", "chromStart", "chromEnd", "name"]) + "\n")
	for interval in binding_regions:
		outf.write("\t".join([chrom, str(interval[0]), str(interval[1]), tf_name]) + "\n")
	outf.close()

	# translate regions if assembly not hg19 or mm9
	if input_assembly != ref_assembly:
		# this is CrossMap now
		outfile_translated = "%s/%s_%s_tfbs_%s.tsv"%(outdir, gene_name if gene_name else chrom, tf_name, ref_assembly)
		(map_tree, target_chrom_sizes, source_chrom_sizes) = read_chain_file(chain_file, print_table = False)
		crossmap_bed_file(map_tree, outfile, outfile_translated)

	return True


#########################################
def main():

	if len(sys.argv) < 2:
		print  ("usage: %s <input_data.tsv>" % sys.argv[0])
		exit()

	outdir  = "raw_data/tf_binding_sites"
	geodir  = "/storage/databases/geo"
	tadfile = "/storage/databases/tads/encode/ENCFF633ORE.bed"
	input_data_file = sys.argv[1]
	for dependency in [outdir, geodir, tadfile, input_data_file]:
		if not os.path.exists(dependency):
			print(dependency,"not found")
			exit()

	inf = open(input_data_file, "r")
	for line in inf:
		if line[0]=='%': continue
		process_tf_binding(outdir, geodir, tadfile, line)

#########################################
########################################
if __name__ == '__main__':
	main()


