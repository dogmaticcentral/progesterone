#!/usr/bin/python3

#
# This file is part of Progesternoe pipeline.
#
# Progesterone pipeline  is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Progesterone pipeline is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Progesterone pipeline.  If not, see <https://www.gnu.org/licenses/>.
#

# the main problem with the local (i.e. downloaded from GEO)
# is the arbitrariness of the format

import sys, os

# pycharm recognizes this if it says .linkto
# however python3 does not like the dot
from utils.mysqldb import *

from utils.utils import *
from utils.CrossMap import *


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

	ucsc_gene_regions_dir = "/storage/databases/ucsc/gene_ranges/%s/%s" % (species, input_assembly)

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
		if species=="human" and input_assembly==ref_assembly: # TODO: fix this
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
	outfile = "%s/%s_%s_tfbs_%s_%s.tsv"%(outdir, gene_name if gene_name else chrom, tf_name, input_assembly, geo_id)
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

	if len(sys.argv) < 3:
		print("Usage: %s <data directory path> <input_data.tsv>" % sys.argv[0])
		print("<data directory path> should be something like \"/storage/databases/geo\"")
		print("and contain bed files with sub-path <experiment id>/<file id>.bed")
		print("<input_data.tsv> should contain tab separated  columns of the form")
		print("organism | gene name | TF name | experiment id | agonist file bed | control/vehicle file bed | assembly")
		print("Control file is optional.")
		exit()

	[datadir, input_data_file] = sys.argv[1:3]
	datadir = datadir.rstrip("/")
	outdir  = "raw_data/tf_binding_sites_%s" % (datadir.split("/").pop())
	tadfile = "/storage/databases/encode/ENCSR551IPY/ENCFF633ORE.bed"
	for dependency in [outdir, datadir, tadfile, input_data_file]:
		if not os.path.exists(dependency):
			print(dependency,"not found")
			exit()

	inf = open(input_data_file, "r")
	for line in inf:
		if line[0]=='%': continue
		process_tf_binding(outdir, datadir, tadfile, line)

#########################################
########################################
if __name__ == '__main__':
	main()


