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

# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
# -A skips auto rehash

from utils.mysqldb import *
import os


#########################################
def main():


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
	species = {'hg18':'human','hg19':'human', 'mm9':'mouse'}

	# note you should have the skip-auto-rehash option in .ucsc_myql_conf
	# it is the equivalent to -A on the mysql command line
	# means: no autocompletion, which makes mysql get up mych faster
	ucsc_db     = connect_to_mysql(ucsc_conf_file)
	ucsc_cursor = ucsc_db.cursor()

	for assembly in ["hg18","hg19", "mm9"]:
		qry = "select chrom, size from %s.chromInfo" % assembly
		rows = search_db(ucsc_cursor,qry)
		if not rows or 'Error' in rows[0][0]:
			search_db(ucsc_cursor,qry, verbose=True)
			break
		for row in rows:
			[chrom, size] = row
			if '_' in chrom: continue  # we don't want to get _too_ general here
			fixed_fields  = {'species':species[assembly], 'chromosome':chrom, 'assembly':assembly, 'rtype':'chromosome'}
			update_fields = {'rfrom':1, 'rto':int(size), 'xref_id':xref_id}
			store_or_update(local_cursor, 'regions', fixed_fields, update_fields)
	ucsc_cursor.close()
	ucsc_db.close()

	local_cursor.close()
	local_db.close()
	return True


#########################################
if __name__ == '__main__':
	main()


