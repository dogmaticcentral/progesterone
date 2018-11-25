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

import math

from utils.mysqldb import *
from utils.utils import *
from statistics import mean, stdev, median
from utils.CrossMap import *


#########################################
def  get_binding_regions_from_ucsc(assembly, chrom, tf_name):
	db     = connect_to_mysql("/home/ivana/.ucsc_mysql_conf")
	cursor = db.cursor()
	switch_to_db(cursor, assembly) # human build name

	table = 'wgEncodeRegTfbsClusteredV3'
	# python thinks these are all strings
	qry = "select chromStart, chromEnd from %s where chrom='chr%s' and name='%s' " % (table, chrom, tf_name)
	# columns: bin, chrom, chromStart, chromEnd, name, score, expCOunt, expNums, expSCores
	ret = search_db(cursor, qry)
	cursor.close()
	db.close()

	if ret==None:
		print ("No ret for %s", qry)
		exit()
	if isinstance([0][0],str) and 'Error'in ret[0][0]:
		print(ret)
		exit()

	return ret


#########################################
def get_binding_regions_from_local_file(ref_assembly, species, chrom, tf_name, input_data_path, input_data_table, scratchdir):
	inf =  open(input_data_table,"r")
	[geo_id, agonist_file, control_file, input_assembly] = [None]*4
	for line in inf:
		[spec, gene_name, tfnm, gid, ag_file, ctrl_file, assm] = line.rstrip().split("\t")
		if spec==species and tfnm==tf_name:
			[geo_id, agonist_file, control_file, input_assembly] = [gid, ag_file, ctrl_file, assm]
			break
	if not geo_id:
		print("Input file not found for", species, tf_name)
		exit()
	binding_regions = read_binding_intervals("{}/{}".format(input_data_path, geo_id), agonist_file, control_file, chrom, None, None)

	if input_assembly == ref_assembly:  # we'll need  tools to translate
		return binding_regions
	else:
		chain_file="/storage/databases/liftover/{}To{}.over.chain".format(input_assembly, ref_assembly.capitalize())
		for dep in [chain_file, scratchdir]:
			if not dep or not os.path.exists(dep):
				print ("I need scratch dir and chain file to translate %s to %s" %(input_assembly, ref_assembly))
				exit()
		outfile = "{}/{}.{}.bed".format(scratchdir, os.getpid(), input_assembly)
		outfile_translated = "{}/{}.{}.bed".format(scratchdir, os.getpid(), ref_assembly)
		outf = open (outfile,"w")
		for interval in binding_regions:
			outf.write("\t".join(["chr%s"%chrom, str(interval[0]), str(interval[1])]) + "\n")
		outf.close()
		(map_tree, target_chrom_sizes, source_chrom_sizes) = read_chain_file(chain_file, print_table = False)
		crossmap_bed_file(map_tree, outfile, outfile_translated)
		inf = open(outfile_translated,"r")
		translated_binding_regions = []
		for line in inf:
			[chrom, start, end] = line.rstrip().split("\t")
			translated_binding_regions.append([start, end])
		os.remove(outfile)
		os.remove(outfile_translated)
		unmapped = outfile_translated+".unmap"
		if os.path.exists(unmapped): os.remove(unmapped)
		return translated_binding_regions


#########################################
def main():

	if len(sys.argv) < 3:
		print  ("usage: %s <species> <tf_name> <chrom> [<data directory path>  <input_data.tsv> <scratch dir>]" % sys.argv[0])
		exit()
	#  this should come from the previous script, 02_emve_tads.py
	[species, tf_name, chrom] = sys.argv[1:4]
	input_data_path  = sys.argv[4] if len(sys.argv) > 4 else None
	input_data_table = sys.argv[5] if len(sys.argv) > 5 else None
	scratchdir = sys.argv[6] if len(sys.argv) > 6 else None

	##############################################
	if species!='human':
		print('as of this writing we only have TADs for human')
		exit()
	##############################################

	ref_assembly = "hg19" if species=="human" else "mm9" # now this is what I call general

	use_local_file =  (species=='mouse' or tf_name=="PGR")
	if use_local_file and (input_data_path==None or input_data_table==None or scratchdir==None):
		print ("to process local  data file I need <data directory path>, <input_data.tsv>, and <scratch dir>")
		exit ()

	outdir = "raw_data"
	tadfile = "/storage/databases/encode/ENCSR551IPY/ENCFF633ORE.bed"
	dependencies = [outdir, tadfile]
	if use_local_file: dependencies += [input_data_path, input_data_table, scratchdir]
	for dep in dependencies:
		if os.path.exists(dep): continue
		print(dep, "not found")
		exit()

	##############################################
	if use_local_file: # we know we don't have this in UCSC
		binding_regions = get_binding_regions_from_local_file(ref_assembly, species, chrom, tf_name,
							input_data_path, input_data_table, scratchdir)
	else:
		binding_regions = get_binding_regions_from_ucsc(ref_assembly, chrom, tf_name)


	##############################################
	tads = get_all_tads (tadfile, "chr%s"%chrom)
	number_of_tads = len(tads)
	bin_counts = [0]*number_of_tads
	for r in binding_regions:
		for t in range(number_of_tads):
			try:
				interval_start = int(r[0])
			except:
				continue
			if tads[t][0]<=interval_start<=tads[t][1]:
				bin_counts[t] += 1
				break
	total_tad_length = 0
	for t in range(number_of_tads):
		total_tad_length += tads[t][1] - tads[t][0] + 1

	##############################################
	avg  = mean(bin_counts)
	stdv = stdev(bin_counts)
	sum_pop = sum(bin_counts)

	outfile = ("{}/{}_{}_{}_tad_distribution.tsv".format(outdir, tf_name, chrom, ref_assembly))
	outf = open(outfile,"w")
	outf.write("\t".join(["% tad_nr", "region", "TFbs count", "TAD length in Kbp",
						"TFbs count fraction", "TFbs count z-score", "TFbs count fraction/TAD length fraction"])+"\n")
	for t in range(number_of_tads):
		z = (bin_counts[t]-avg)/stdv
		bin_fraction =  bin_counts[t]/sum_pop
		tad_fraction = (tads[t][1] - tads[t][0] + 1)/total_tad_length
		outf.write("%d \t%s\t%d\t%d\t%.3f\t%.2f\t%.2f\n" %
					(t,  "chr{}:{}-{}".format(chrom, tads[t][0], tads[t][1]),
					bin_counts[t],  int((tads[t][1] - tads[t][0] + 1)/1000),
					bin_fraction, z, bin_fraction/tad_fraction) )
	outf.close()


	return True


#########################################
########################################
if __name__ == '__main__':
	main()

