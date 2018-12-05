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

	local_conf_file = "/home/ivana/.mysql_conf"
	ucsc_conf_file  = "/home/ivana/.ucsc_mysql_conf"

	for dependency in [local_conf_file, ucsc_conf_file]:
		if not os.path.exists(dependency):
			print(dependency, "not found")
			exit()

	ucsc_db     = connect_to_mysql(ucsc_conf_file)
	ucsc_cursor = ucsc_db.cursor()
	switch_to_db(ucsc_cursor, 'hgcentral')
	qry = "select organism, nibPath, scientificName from dbDb "
	ucsc_ret = search_db(ucsc_cursor,qry)
	hard_check(ucsc_db, ucsc_cursor, ucsc_ret, qry)
	ucsc_cursor.close()
	ucsc_db.close()


	local_db = connect_to_mysql(local_conf_file)
	local_cursor = local_db.cursor()
	# autocommit is on by default, except when it is not
	search_db(local_cursor,"set autocommit=1")
	switch_to_db(local_cursor,'progesterone')
	for [organism, path, sciname] in ucsc_ret:
		assm = path.split("/")[2]
		organism = organism.replace("'","")
		print(organism, assm, sciname)
		fields = {'common_name':organism,'scientific_name':sciname,  'assembly':assm}
		store_or_update(local_cursor,'assemblies',fields, None)

	local_cursor.close()
	local_db.close()
	return True


#########################################
if __name__ == '__main__':
	main()


