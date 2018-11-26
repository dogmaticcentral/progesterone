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

from utils.mysqldb import *
import os


#########################################
def main():

	species ="mouse"
	assembly = "mm9"
	storage = "/storage/databases/ucsc/gene_ranges/%s/%s" % (species,assembly)
	if not os.path.exists(storage):
		print("Please create %s" % storage)
		exit()

	db     = connect_to_mysql("/home/ivana/.ucsc_mysql_conf")
	cursor = db.cursor()

	switch_to_db(cursor, assembly)
	chromosomes = []
	if species =='human':
		chromosomes = ["chr"+str(x) for x in range(1,23)] + ["chrX"]
	elif species =='mouse':
		chromosomes = ["chr"+str(x) for x in range(1,20)] + ["chrX"]
	for chrom in chromosomes:
		print("downloading data for", chrom)
		outf = open("/storage/databases/ucsc/gene_ranges/{}/{}/{}.csv".format(species,assembly,chrom), "w")
		outf.write( "\t".join(["name", "name2", "strand", "txStart", "txEnd"]) )
		qry  = "select name,  name2, strand, txStart, txEnd "
		qry += "from refGene "
		qry += "where chrom='%s' " % chrom
		qry += "and name like 'NM_%'"   # refseq says: NM_	mRNA	Protein-coding transcripts (usually curated)
		rows = search_db(cursor,qry)
		for row in rows:
			outf.write("\t".join( [ str(r) for r in row])+"\n")
		outf.close()
	cursor.close()
	db.close()
	return True


#########################################
if __name__ == '__main__':
	main()


