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

	gene_name = "Hand2"
	assembly = "hg19"
	external_exp_id = "ENCFF633ORE"
	conf_file = "/home/ivana/.mysql_conf"

	for prerequisite in [ conf_file]:
		if os.path.exists(prerequisite): continue
		print(prerequisite, "not found")
		exit()

	db = connect_to_mysql(conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,'progesterone')

	# find xref_id for the experimental data file
	exp_file_xref_id = get_xref_id(db,cursor,external_exp_id)

	# find gene coordinates
	[chromosome, strand, min_start, max_end] = get_gene_coords(db,cursor,gene_name,assembly)

	# finally, use that info to find the TAD
	[tad_start, tad_end] = get_tad_region(db, cursor, exp_file_xref_id, chromosome, min_start, max_end)

	cursor.close()
	db.close()

	print ( "{} {} {}:{}-{}".format(gene_name, strand, chromosome, min_start, max_end) )
	print ("TAD containing %s region: %s:%d-%d   length %d"%(gene_name, chromosome,
															 tad_start, tad_end, tad_end-tad_start+1))

#########################################
if __name__ == '__main__':
	main()




