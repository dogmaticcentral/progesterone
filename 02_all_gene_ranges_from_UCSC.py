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

	species ="human"
	assembly = "hg19"

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

	ucsc_db     = connect_to_mysql(ucsc_conf_file)
	ucsc_cursor = ucsc_db.cursor()
	switch_to_db(ucsc_cursor, assembly)

	chromosomes = []
	if species =='human':
		chromosomes = ["chr"+str(x) for x in range(1,23)] + ["chrX", "chrY"]
	elif species =='mouse':
		chromosomes = ["chr"+str(x) for x in range(1,20)] + ["chrX", "chrY"]
	for chrom in chromosomes:
		print("downloading data for", assembly, chrom)
		qry  = "select name,  name2, strand, txStart, txEnd "
		qry += "from refGene "
		qry += "where chrom='%s' " % chrom
		qry += "and name like 'NM_%'"   # refseq says: NM_	mRNA	Protein-coding transcripts (usually curated)
		rows = search_db(ucsc_cursor,qry)
		print("loading ...")
		for row in rows:
			[gene_name, strand, rfrom, rto] = row[1:]
			# store region
			fields  = {'species':species, 'chromosome':chrom, 'assembly':assembly, 'rtype':'gene',
						'rfrom':rfrom, 'rto':rto, 'strand':strand, 'xref_id':xref_id}
			region_id = store_without_checking(local_cursor, 'regions', fields)
			if region_id<0:
				print("insert failure for", gene_name, strand, rfrom, rto)
				exit()
			# store gene
			fields  = {'name':gene_name, 'region_id':region_id}
			gene_id = store_without_checking(local_cursor, 'genes', fields)

	ucsc_cursor.close()
	ucsc_db.close()


	local_cursor.close()
	local_db.close()
	return True


#########################################
if __name__ == '__main__':
	main()


