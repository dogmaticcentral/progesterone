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
def transform_coords(scratchdir, chain_file, chrom, binding_regions):

	outfile = "%s/tfbs_%s.tsv"%(scratchdir, os.getpid())
	outf = open (outfile,"w")
	for interval in binding_regions:
		outf.write("\t".join([chrom, str(interval[0]), str(interval[1])]) + "\n")
	outf.close()

	# this is CrossMap now
	outfile_translated  = "%s/tfbs_translated_%s.tsv"%(scratchdir, os.getpid())
	(map_tree, target_chrom_sizes, source_chrom_sizes) = read_chain_file(chain_file, print_table = False)
	crossmap_bed_file(map_tree, outfile, outfile_translated)

	#read binding regions back in
	with open(outfile_translated,"r") as inf:
		new_binding_regions = [line.split("\t")[1:3] for line in inf.read().split("\n") if len(line.replace(" ",""))>0]
	# remove aux files
	os.remove(outfile)
	os.remove(outfile_translated)

	return new_binding_regions

#########################################
def get_binding_regions (db, cursor, scratchdir, data_homedir, input_line):

	#############################
	# input processing
	[species, gene_or_chr, tad_exp_id, tf_name, source, datafile_id, agonist_file, ctrl_file, input_assembly] = input_line.rstrip().split("\t")

	chain_file_input_to_ref_assmb =None
	chain_file_ref_to_input_assmb =None

	ref_assembly = None
	if species=="human":
		ref_assembly="hg19"
	elif species=="mouse":
		ref_assembly="mm9"
	else:
		print("organism \"%s\" not recognized"%species)
		exit()

	if input_assembly != ref_assembly:  # we'll need  tools to translate
		chain_file_input_to_ref_assmb ="/storage/databases/liftover/{}To{}.over.chain".format(input_assembly, ref_assembly.capitalize())
		chain_file_ref_to_input_assmb ="/storage/databases/liftover/{}To{}.over.chain".format(ref_assembly, input_assembly.capitalize())

	if ctrl_file.replace(" ","") == "": ctrl_file = None
	data_dir = "{}/{}".format(data_homedir, datafile_id)

	dependencies = [data_dir, "{}/{}".format(data_dir,agonist_file)]
	if ctrl_file: dependencies += ["{}/{}".format(data_dir,ctrl_file)]
	if chain_file_input_to_ref_assmb : dependencies += [chain_file_input_to_ref_assmb ]
	for dep in dependencies:
		if not os.path.exists(dep):
			print(dep,"not found")
			exit()

	if 'chr' in gene_or_chr:
		chrom = gene_or_chr
		region_start =  region_end = None
	else:
		gene_name = gene_or_chr
		[chrom, strand, min_start, max_end] = get_gene_coords(db,cursor, gene_name, ref_assembly)
		if tad_exp_id==None or tad_exp_id.replace(" ","") == "":
			region_start = min_start - 1.0e6
			region_end   = max_end   + 1.0e6
		else:
			exp_file_xref_id = get_xref_id(db,cursor,tad_exp_id)
			[region_start, region_end] = get_tad_region(db, cursor, exp_file_xref_id, chrom, min_start, max_end)
		if input_assembly != ref_assembly:
			[[region_start, region_end]] = transform_coords(scratchdir, chain_file_ref_to_input_assmb ,
															chrom, [[region_start, region_end]])

	#############################
	# business
	binding_regions = read_binding_intervals(data_dir, agonist_file, ctrl_file, chrom.replace('chr', ''), region_start, region_end)
	#############################
	# transform coords if needed
	if input_assembly != ref_assembly:
		binding_regions = transform_coords(scratchdir, chain_file_input_to_ref_assmb , chrom, binding_regions)

	return [tf_name, source, datafile_id, species, chrom, ref_assembly, binding_regions]

#########################################
def store_binding_regions(cursor, binding_info, check_duplicates = False):
	[tf_name, source, datafile_id, species, chrom, assembly, binding_regions] = binding_info
	xref_id = store_xref(cursor, source, datafile_id)
	for [start,end] in binding_regions:
		# store region (address)
		fields  = {'species':species, 'chromosome':chrom, 'assembly':assembly, 'rtype':'chipseq',
					'rfrom':start, 'rto':end, 'xref_id':xref_id}
		region_id = store_or_update(cursor, 'regions', fields, None) if check_duplicates else \
					store_without_checking(cursor, 'regions', fields)

		# store info about the binding region
		fields = {'tf_name':tf_name, 'chipseq_region_id':region_id, 'xref_id':xref_id}
		binding_site_id = store_or_update (cursor, 'binding_sites', fields, None) if check_duplicates else \
						store_without_checking(cursor, 'binding_sites', fields)
	return

#########################################
def main():

	check_duplicates = True # when storing to the local database - makes things slower

	if len(sys.argv) < 3:
		print("Usage: %s <data directory path> <input_data.tsv> " % sys.argv[0])
		print("<data directory path> should be something like \"/storage/databases/geo\"")
		print("and contain bed files with sub-path <experiment id>/<file id>.bed")
		print("<input_data.tsv> should contain tab separated  columns of the form")
		print("organism |  gene or chr name | TAD exp id | TF name | source | experiment id | agonist file bed | control/vehicle file bed | assembly")
		print("Source refers to the source database. Control file is optional, and so is TAD experiment ID.")
		print("If control file is given, peaks will be subtracted from agonist peaks.")
		print("TAD exp is ignored if a chromosome name is given. If TAD exp id is given,")
		print("TAD region encompassing the gene will be used as target range.")
		print("Otherwise 1Mbp on each side of the gene will be used. ")
		exit()

	[data_homedir, input_data_file] = sys.argv[1:3]
	data_homedir = data_homedir.rstrip("/")
	scratchdir   = "/home/ivana/scratch"
	local_conf_file = "/home/ivana/.mysql_conf"
	for dependency in [data_homedir, input_data_file, local_conf_file, scratchdir]:
		if not os.path.exists(dependency):
			print(dependency,"not found")
			exit()

	db = connect_to_mysql(local_conf_file)
	cursor = db.cursor()
	search_db(cursor,"set autocommit=1")
	switch_to_db(cursor,'progesterone')

	inf = open(input_data_file, "r")
	for line in inf:
		if line[0]=='%': continue
		if len(line.replace(" ","").replace("\t",""))==0: continue
		line = line.rstrip()
		binding_info = get_binding_regions(db, cursor, scratchdir, data_homedir, line)
		store_binding_regions(cursor, binding_info, check_duplicates)
	cursor.close()
	db.close()

#########################################
########################################
if __name__ == '__main__':
	main()


