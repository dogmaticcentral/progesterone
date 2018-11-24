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

# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
# -A skips auto rehash
from utils.mysqldb import *
import os


#########################################
def main():

	storage = "/storage/databases/ucsc/chromosome_lengths"
	if not os.path.exists(storage):
		print("Please create %s." % storage)
		exit()

	# note the skip-auto-rehash option in .ucsc_myql_conf
	# it is the equivalent to -A on the mysql command line
	# means: no autocompletion, which makes mysql get up mych faster
	db     = connect_to_mysql("/home/ivana/.ucsc_mysql_conf")
	cursor = db.cursor()

	for assembly in ["hg18","hg19", "mm9"]:
		qry = "select chrom, size from %s.chromInfo" % assembly
		rows = search_db(cursor,qry)
		if not rows or 'Error' in rows[0][0]:
			search_db(cursor,qry, verbose=True)
			break
		outf = open("%s/%s.tsv" % (storage,assembly),"w")
		for row in rows:
			outf.write("\t".join( [ str(r) for r in row])+"\n")
		outf.close()

	cursor.close()
	db.close()
	return True


#########################################
if __name__ == '__main__':
	main()


