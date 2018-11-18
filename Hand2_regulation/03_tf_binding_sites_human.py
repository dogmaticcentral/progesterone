#!/usr/bin/python3
# mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A
# -A skips auto rehash
import sys

# pycharm recognizes this if it says .linkto
# however python3 does not like the dot
from linkto_python_modules.mysqldb import *

# finding exps by exp num:
#  SET @row_number = 0;  select * from (SELECT (@row_number:=@row_number + 1) \
#  AS expNum, source, factor, cellType, treatment, lab from wgEncodeRegTfbsClusteredInputsV3) as t where factor='ESR1';
# (the wgEncodeRegTfbsClusteredInputsV3 table has no id nums, so we have to count the rows)
#########################################
def main():

	assembly = "hg19"

	if len(sys.argv) < 5:
		print  ("usage: %s <gene_name> <chrom>  <from>  <to>" % sys.argv[0])
		exit()
	#  this should come from the previous script, 02_emve_tads.py
	[gene_name, chrom, start,end] = sys.argv[1:5]
	db     = connect_to_mysql("/home/ivana/.ucsc_mysql_conf")
	cursor = db.cursor()
	switch_to_db(cursor, assembly) # human build name

	# our table du jour is wgEncodeRegTfbsClusteredV3;
	table = 'wgEncodeRegTfbsClusteredV3'
	# python thinks these are all strings
	qry = "select * from %s where chrom='chr%s' and chromStart>%s and  chromEnd<%s" % (table, chrom, start,end)
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

	outf = open ("raw_data/%s_tfbs_%s.tsv"%(gene_name,assembly),"w")
	outf.write("\t".join(["% chrom", "chromStart", "chromEnd", "name", "score", "expCount", "expNums", "expSCores"]) + "\n")
	outf.write("\n".join( "\t".join([str(field).replace("b'","").replace("'","") for field in row[1:]]) for row in ret) + "\n")
	outf.close()

	return True


#########################################
########################################
if __name__ == '__main__':
	main()

