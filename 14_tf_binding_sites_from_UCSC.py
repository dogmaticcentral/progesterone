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

from utils.mysqldb import  *
import os, sys


#########################################
def main():

	assembly = "hg19" # afaik this is the only assembly with ENCODE data
	outdir = "raw_data/tf_binding_sites"
	if not os.path.exists(outdir):
		print(outdir,"not found (make that  dir for output, or change %s appropriately)" % sys.argv[0])
		exit()

	#  the broader chromosome region (such as TAD)
	#  should come from the previous script, 02_emve_tads.py
	if len(sys.argv) < 5:
		print  ("usage:   %s <gene_name> <chrom>  <from>  <to>" % sys.argv[0])
		print  ("example: %s  Hand2 4 173880001 175320000" % sys.argv[0])
		exit()

	[gene_name, chrom, start,end] = sys.argv[1:5]
	db     = connect_to_mysql("/home/ivana/.ucsc_mysql_conf")
	cursor = db.cursor()
	switch_to_db(cursor, assembly) # human build name

	# our table du jour is wgEncodeRegTfbsClusteredV3;
	table = 'wgEncodeRegTfbsClusteredV3'
	# python thinks these are all strings
	qry = "select * from %s where chrom='chr%s' and chromStart>%s and  chromEnd<%s" % (table, chrom, start,end)
	# columns: bin, chrom, chromStart, chromEnd, name, score, expCount, expNums, expScores
	ret = search_db(cursor, qry)

	cursor.close()
	db.close()


	if ret==None:
		print ("No ret for %s", qry)
		exit()
	if isinstance([0][0],str) and 'Error'in ret[0][0]:
		print(ret)
		exit()

	outf = open ("%s/%s_tfbs_%s.tsv"%(outdir, gene_name,assembly),"w")
	outf.write("\t".join(["% chrom", "chromStart", "chromEnd", "name", "score", "expCount", "expNums", "expSCores"]) + "\n")
	outf.write("\n".join( "\t".join([str(field).replace("b'","").replace("'","") for field in row[1:]]) for row in ret) + "\n")
	outf.close()

	return True


#########################################
########################################
if __name__ == '__main__':
	main()

