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

from utils.utils import *
from utils.mysqldb import *


#########################################
def main():

	species = ['human','mouse']
	tf_names = ['PGR', 'ESR1']

	conf_file  = "/home/ivana/.mysql_conf"

	db = connect_to_mysql(conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,'progesterone')

	for sp in species:
		print()
		print(sp)
		for tf_name in tf_names:
			print("\t",tf_name)

			# number of binding sites
			qry = "select count(1) from binding_sites as b, regions as r "
			qry += "where b.tf_name='%s' " % tf_name
			qry += "and r.species='%s' " % sp
			qry += "and b.region_id=r.id "
			number_of_bs = search_db(cursor,qry)[0][0]
			print("number of (chipseq) binding sites:", number_of_bs)

			# number of binding sites with motif assigned
			subqry = "select count(1) from  binding_site2motif where binding_site_id=b.id"
			qry += "and (%s)>0 " % subqry
			bs_w_motif = search_db(cursor,qry)[0][0]
			print("number of binding sites w motif:", bs_w_motif)

			# number of binding sites with motif that has alignment associated
			qry  = "select count(1) from alignments as a, motifs as m, regions as r "
			qry += "where a.id=m.alignment_id and m.region_id=r.id "
			qry += "and m.tf_name='%s' and r.species='%s'" % (tf_name, sp)
			motif_w_almt = search_db(cursor,qry)[0][0]
			print("number of  motifs w alignment:", motif_w_almt)


	cursor.close()
	db.close()


#########################################
########################################
if __name__ == '__main__':
	main()
