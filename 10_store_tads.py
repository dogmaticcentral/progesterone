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


# the single input file from endometrial microvascular endothelial cells from
# Job Dekker lab, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105710
# the (bed) file with TADs can be found here
# https://www.encodeproject.org/experiments/ENCSR551IPY/
# (under processed data)

from utils.utils import *
from utils.mysqldb import *

#########################################
def main():

	species   = "human"
	assembly  = "hg19" # this came from the file metadata, no the file itself
	tadfile   = "/storage/databases/encode/ENCSR551IPY/ENCFF633ORE.bed"
	conf_file = "/home/ivana/.mysql_conf"
	########
	# references
	# I don't see  a way to automate this - for example we might have had the pubmed ref, but we don't
	geo_exp_id  = 'ENCSR551IPY'
	geo_file_id = 'ENCFF633ORE'

	for prerequisite in [tadfile, conf_file]:
		if os.path.exists(prerequisite): continue
		print(prerequisite, "not found")
		exit()

	tads = {}
	inf = open(tadfile, "r")
	for line in inf:
		[chr, start, end] = line.rstrip().split()[:3]
		if not chr in tads: tads[chr] = []
		tads[chr].append([int(start), int(end)])

	db = connect_to_mysql(conf_file)
	cursor = db.cursor()
	search_db(cursor,"set autocommit=1")
	switch_to_db(cursor,'progesterone')
	xref_exp_id  = store_xref(cursor, 'geo', geo_exp_id)
	xref_file_id = store_xref(cursor, 'geo', geo_file_id, parent_id=xref_exp_id)

	for chr, regions in tads.items():
		for [start, end] in regions:
			# store region
			fields  = {'species':species, 'chromosome':chr, 'assembly':assembly, 'rtype':'tad',
						'rfrom':start, 'rto':end,  'xref_id':xref_file_id}
			region_id = store_without_checking(cursor, 'regions', fields)

	cursor.close()
	db.close()


#########################################
if __name__ == '__main__':
	main()




