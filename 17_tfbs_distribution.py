#!/usr/bin/python3

#
# This file is part of Progesterone pipeline.
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
def main():

	if len(sys.argv) < 4:
		print  ("usage: %s <species> <tf_name> <chrom>" % sys.argv[0])
		exit()
	[species, tf_name, chromosome] = sys.argv[1:4]
	if not 'chr' in chromosome: chromosome = 'chr'+chromosome

	##############################################
	if species!='human':
		print('as of this writing we only have TADs for human')
		exit()
	##############################################

	ref_assembly = "hg19" if species=="human" else "mm9" # now this is what I call general

	external_exp_id = "ENCFF633ORE" # the only HiC experiment defining TADs we are using so far

	conf_file = "/home/ivana/.mysql_conf"
	outdir = "raw_data/tad_distributions"
	dependencies = [outdir, conf_file]
	for dep in dependencies:
		if os.path.exists(dep): continue
		print(dep, "not found")
		exit()

	##############################################
	db = connect_to_mysql(conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,'progesterone')
	exp_file_xref_id = get_xref_id(db,cursor,external_exp_id)

	binding_regions = get_binding_regions(db, cursor, ref_assembly, chromosome, tf_name)
	tads = get_all_tads (db, cursor, exp_file_xref_id, chromosome)
	cursor.close()
	db.close()



	##############################################
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

	outfile = ("{}/{}_{}_{}_tad_distribution.tsv".format(outdir, tf_name, chromosome, ref_assembly))
	outf = open(outfile,"w")
	outf.write("\t".join(["% tad_nr", "region", "TFbs count", "TAD length in Kbp",
						"TFbs count fraction", "TFbs count z-score", "TFbs count fraction/TAD length fraction"])+"\n")
	for t in range(number_of_tads):
		z = (bin_counts[t]-avg)/stdv
		bin_fraction =  bin_counts[t]/sum_pop
		tad_fraction = (tads[t][1] - tads[t][0] + 1)/total_tad_length
		outf.write("%d\t%s\t%d\t%d\t%.3f\t%.2f\t%.2f\n" %
					(t,  "{}:{}-{}".format(chromosome, tads[t][0], tads[t][1]),
					bin_counts[t],  int((tads[t][1] - tads[t][0] + 1)/1000),
					bin_fraction, z, bin_fraction/tad_fraction) )
	outf.close()


	return True


#########################################
########################################
if __name__ == '__main__':
	main()

