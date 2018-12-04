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
# UCSC doe not have the pointers directly back to ENCODE, so I found them for ESR1 - by hand
# encode_esr1_xps.tsv must contain 3 columns: UCSC id, encode experiment id, and encode file id
#########################################
def main():

	conf_file    = "/home/ivana/.mysql_conf"
	mapping_file = "encode_esr1_xps.tsv"
	for dependency in [conf_file, mapping_file]:
		if not os.path.exists(dependency):
			print(dependency,"not found")
			exit()
	encode_exp_id  = {}
	encode_file_id = {}
	ucsc_ids = []
	with open(mapping_file,"r") as inf:
		for line in inf:
			if 'UCSC' in line: continue # header
			[ucsc, encode_exp, encode_file] = line.split("\t")[:3]
			ucsc_ids.append(ucsc)
			encode_exp_id[ucsc] = encode_exp
			encode_file_id[ucsc]   = encode_file

	#########################
	# plug in to local database
	db = connect_to_mysql(conf_file)
	cursor = db.cursor()
	search_db(cursor,"set autocommit=1")
	switch_to_db(cursor,'progesterone')

	# this might not be the best idea if the database grows really large
	# first make sure we have single entry for each of multiple ids
	for line in search_db(cursor,"select id, external_id from xrefs where xtype='ucsc'"):
		[xref_id, ucsc_str] = line
		ucsc_ids_stored = ucsc_str.split(",")
		if len(ucsc_ids_stored) <2: continue
		for ucsc_id in ucsc_ids_stored:
			store_or_update(cursor, 'xrefs', {'xtype':'ucsc', 'external_id':ucsc_id}, None)

	# now for each single entry, make parent point to encode file, and encode file's parent to encode exp
	for line in search_db(cursor,"select id, external_id from xrefs where xtype='ucsc' and external_id not like '%,%'"):
		[ucsc_xref_id, ucsc_id] = line
		if not ucsc_id in ucsc_ids: continue
		encode_file_xref_id = store_or_update(cursor, 'xrefs', {'xtype':'encode', 'external_id': encode_file_id[ucsc_id]}, None)
		search_db(cursor, "update xrefs set parent_id=%d where id=%d" % (encode_file_xref_id, ucsc_xref_id))
		encode_exp_xref_id = store_or_update(cursor, 'xrefs', {'xtype':'encode', 'external_id': encode_exp_id[ucsc_id]}, None)
		search_db(cursor, "update xrefs set parent_id=%d where id=%d" % (encode_exp_xref_id, encode_file_xref_id))

	cursor.close()
	db.close()


	return True


#########################################
########################################
if __name__ == '__main__':
	main()

