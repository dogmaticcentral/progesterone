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

from utils.mysqldb import *
import os

#########################################
def main():

	gene_name = "Hand2"
	species = ['human', 'gorilla','rhesus',  'mouse', 'rat']

	# The UCSC Genome Browser database: 2019 update.: pubmed id 30407534
	pubmed_id = '30407534'
	local_conf_file = "/home/ivana/.mysql_conf"
	ucsc_conf_file  = "/home/ivana/.ucsc_mysql_conf"

	for dependency in [local_conf_file, ucsc_conf_file]:
		if not os.path.exists(dependency):
			print(dependency, "not found")
			exit()

	local_db = connect_to_mysql(local_conf_file)
	local_cursor = local_db.cursor()
	# autocommit is on by default, except when it is not
	search_db(local_cursor,"set autocommit=1")
	switch_to_db(local_cursor,'progesterone')
	# store reference info
	xref_id = store_xref(local_cursor, 'pubmed', pubmed_id)
	# find all assemblies corresponding to our sample species
	assemblies = {}
	for sp in species:
		qry = "select assembly from assemblies where common_name='%s' " % sp
		ret = search_db(local_cursor,qry)
		hard_check(local_db, local_cursor, ret, qry)
		assemblies[sp] = [r[0] for r in ret]

	ucsc_db     = connect_to_mysql(ucsc_conf_file)
	ucsc_cursor = ucsc_db.cursor()

	for sp in species:
		for assembly in assemblies[sp]:
			print("\nlooking for %s in %s" % (gene_name, assembly))
			ret = switch_to_db(ucsc_cursor, assembly)
			if not ret: continue
			qry  = "select  chrom, strand, txStart, txEnd "
			qry += "from refGene "
			qry += "where name2='%s' " % gene_name
			rows = search_db(ucsc_cursor,qry)
			if not rows:
				qry = qry.replace("refGene","xenoRefGene")
				rows = search_db(ucsc_cursor,qry)
				if not rows: continue

			print(rows)
			if not type(rows[0][0])==str or rows[0][0]=='Error': continue
			print("loading ...")
			for row in rows:
				[chrom, strand, rfrom, rto] = row
				# store region
				fields  = {'species':sp, 'chromosome':chrom, 'assembly':assembly, 'rtype':'gene',
							'rfrom':rfrom, 'rto':rto, 'strand':strand, 'xref_id':xref_id}
				region_id = store_or_update(local_cursor, 'regions', fields, None)
				if region_id<0:
					print("insert failure for", gene_name, strand, rfrom, rto)
					exit()
				# store gene
				fields  = {'name':gene_name, 'region_id':region_id}
				gene_id = store_or_update(local_cursor, 'genes', fields, None)

	ucsc_cursor.close()
	ucsc_db.close()


	local_cursor.close()
	local_db.close()
	return True


#########################################
if __name__ == '__main__':
	main()


