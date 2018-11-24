#!/usr/bin/python3
# mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A
# -A skips auto rehash
import sys, math, os

# pycharm recognizes this if it says .linkto
# however python3 does not like the dot
from utils.mysqldb import *
from utils.utils import *
from statistics import mean, stdev, median

#########################################
def arbitrary_binning(ret):
	number_of_bins = 100
	minpos = min([r[0] for r in ret])
	maxpos = max([r[1] for r in ret])
	binning_range= [int(math.floor(minpos/1000))*1000, int(math.ceil(maxpos/1000))*1000]
	binning_size = (binning_range[1]-binning_range[0])/number_of_bins

	bin_counts = [0]*number_of_bins
	for r in ret:
		i  = math.floor((r[0] - binning_range[0])/binning_size)
		bin_counts[i] += 1

	for i in range(number_of_bins):
		print(i, bin_counts[i])

#########################################
def main():

	# human hg19 hardcoded - I do not know TADs for a general case

	if len(sys.argv) < 3:
		print  ("usage: %s <transcr_factor_name> <chrom> " % sys.argv[0])
		exit()
	#  this should come from the previous script, 02_emve_tads.py
	[tf_name, chrom] = sys.argv[1:3]

	tadfile = "/storage/databases/tads/encode/ENCFF633ORE.bed"
	dependencies = [tadfile]
	tfbs_file = None
	if tf_name=="PGR":
		tfbs_file = "raw_data/tf_binding_sites/chr{}_{}_tfbs_hg19.tsv".format(chrom, tf_name)
		dependencies += [tfbs_file]
	for dep in dependencies:
		if os.path.exists(dep): continue
		print(dep, "not found")
		exit()

	if tf_name=="PGR": # we know we don't have this in UCSC
		with open(tfbs_file,"r") as inf:
			ret = [line.split("\t")[1:3] for line in inf.read().split("\n")]
	else:
		db     = connect_to_mysql("/home/ivana/.ucsc_mysql_conf")
		cursor = db.cursor()
		switch_to_db(cursor, "hg19") # human build name

		# our table du jour is wgEncodeRegTfbsClusteredV3;
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


	tads = get_all_tads (tadfile, "chr%s"%chrom)
	number_of_tads = len(tads)
	bin_counts = [0]*number_of_tads
	for r in ret:
		for t in range(number_of_tads):
			try:
				interval_start = int(r[0])
			except:
				continue
			if tads[t][0]<=interval_start<=tads[t][1]:
				bin_counts[t] += 1
				break

	avg  = mean(bin_counts)
	stdv = stdev(bin_counts)
	med = median(bin_counts)
	max_pop = max(bin_counts)
	sum_pop = sum(bin_counts)
	#print(avg, stdv, med)
	#print()
	for t in range(number_of_tads):
		#if  bin_counts[t]> 15:

		z = "%3.1f"%( (bin_counts[t]-avg)/stdv )
		print ("{} chr{}:{}-{}  {} {}".format(t, chrom, tads[t][0], tads[t][1],  bin_counts[t]/sum_pop, z))



	return True


#########################################
########################################
if __name__ == '__main__':
	main()

