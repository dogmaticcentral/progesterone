#!/usr/bin/python3
# mysql --user=genome --host=genome-mysql.soe.ucsc.edu -A
# -A skips auto rehash
import os

# pycharm recognizes this if it says .linkto
# however python3 does not like the dot
from utils.mysqldb import *

# finding exps by exp num:
#  SET @row_number = 0;  select * from (SELECT (@row_number:=@row_number + 1) \
#  AS expNum, source, factor, cellType, treatment, lab from wgEncodeRegTfbsClusteredInputsV3) as t where factor='ESR1';
# (the wgEncodeRegTfbsClusteredInputsV3 table has no id nums, so we have to count the rows)
#########################################
def main():

	assembly = "hg19" # afaik this is the only assembly with ENCODE data
	if len(sys.argv) < 2:
		print  ("usage:   %s <tf name> " % sys.argv[0])
		exit()

	tf_name = sys.argv[1]
	db     = connect_to_mysql("/home/ivana/.ucsc_mysql_conf")
	cursor = db.cursor()
	switch_to_db(cursor, assembly) # human build name

	# our table du jour is wgEncodeRegTfbsClusteredV3;
	table = 'wgEncodeRegTfbsClusteredInputsV3'
	# python thinks these are all strings
	subquery  = "SELECT (@row_number:=@row_number + 1) "
	subquery += "AS expNum, source, factor, cellType, treatment, lab from %s " % table
	qry = "SET @row_number = 0;  select * from (%s) as t  where factor='%s';" % (subquery, tf_name)
	# columns: bin, chrom, chromStart, chromEnd, name, score, expCount, expNums, expScores
	print(qry)
	ret = cursor.executemany( """"SET @row_number = 0;  select * from (SELECT (@row_number:=@row_number + 1)AS expNum, source, factor, cellType, treatment, lab from wgEncodeRegTfbsClusteredInputsV3) as t where factor='ESR1';"""
	                          , [])


	if ret==None:
		print ("No ret for %s" %qry)
		exit()
	if not ret or  isinstance([0][0],str) and 'Error'in ret[0][0]:
		search_db(cursor, qry, verbose=True)
		exit()

	print(ret)

	cursor.close()
	db.close()


	return True


#########################################
########################################
if __name__ == '__main__':
	main()

