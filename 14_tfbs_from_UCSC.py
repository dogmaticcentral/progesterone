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

from utils.mysqldb import  *
import os, sys


#########################################
def main():

	check_duplicates = True # when storing to the local database - makes things slower

	assembly = "hg19" # afaik this is the only assembly with ENCODE data
	species="human"
	ucsc_conf_file  = "/home/ivana/.ucsc_mysql_conf"
	local_conf_file = "/home/ivana/.mysql_conf"

	for prerequisite in [ucsc_conf_file, local_conf_file]:
		if os.path.exists(prerequisite): continue
		print(prerequisite, "not found")
		exit()

	#  input - a  broader chromosome region (such as TAD)
	if len(sys.argv) < 5:
		print  ("usage:   %s  <chrom> <tf_name>  <from|'none'>  <to|'none'>" % sys.argv[0])
		print  ("example: %s  4 ESR1 173880001 175320000" % sys.argv[0])
		exit()

	[chrom, tf_name, start, end] = sys.argv[1:5]

	print ("downloading from ucsc ...")
	db     = connect_to_mysql(ucsc_conf_file)
	cursor = db.cursor()
	switch_to_db(cursor, assembly) # human build name
	# our table du jour is wgEncodeRegTfbsClusteredV3;
	table = 'wgEncodeRegTfbsClusteredV3'
	# python thinks these are all strings
	qry = "select * from %s where chrom='chr%s' and name='%s' " % (table, chrom, tf_name)
	if start!='none': qry += "and chromStart>%s " % start
	if end!='none':   qry += "and chromEnd<%s" % end
	# columns: bin, chrom, chromStart, chromEnd, name, score, expCount, expNums, expScores
	ucsc_ret = search_db(cursor, qry)
	hard_check (db,cursor, ucsc_ret, qry)
	cursor.close()
	db.close()

	print("loading to local db ...")
	db = connect_to_mysql(local_conf_file)
	cursor = db.cursor()
	search_db(cursor,"set autocommit=1")
	switch_to_db(cursor,'progesterone')
	for row in ucsc_ret:

		[bin, chrom, chromStart, chromEnd, name, score, expCount, expNums, expScores] = row
		# store reference
		expid = ",".join([str(i) for i in sorted([int(s) for s in expNums.decode("utf-8").replace(" ","").split(",")])])
		if len(expid)>255: expid='many' # 255 is the storage size for this field - I am not sure what's with all the refs in some cases
		xref_id = store_xref(cursor, 'ucsc', expid)

		# store region (address)
		fields  = {'species':species, 'chromosome':chrom, 'assembly':assembly, 'rtype':'chipseq',
						'rfrom':chromStart, 'rto':chromEnd, 'xref_id':xref_id}
		region_id = store_or_update(cursor, 'regions', fields, None) if check_duplicates else \
					store_without_checking(cursor, 'regions', fields)

		# store info about the binding region
		fields = {'tf_name':name, 'region_id':region_id, 'xref_id':xref_id}
		binding_site_id = store_or_update (cursor, 'binding_sites', fields, None) if check_duplicates else \
						store_without_checking(cursor, 'binding_sites', fields)


	cursor.close()
	db.close()

	return True


#########################################
########################################
if __name__ == '__main__':
	main()

