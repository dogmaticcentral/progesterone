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


# the single input file from endometrial microvascular endothelial cells from
# Job Dekker lab, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105710
# the (bed) file with TADs can be found here
# https://www.encodeproject.org/experiments/ENCSR551IPY/
# (under processed data)

from utils.utils import *
from utils.mysqldb import *

#########################################
def main():

	gene_name = "TP53"
	assembly = "hg19"
	tadfile = "/storage/databases/encode/ENCSR551IPY/ENCFF633ORE.bed"
	conf_file = "/home/ivana/.mysql_conf"

	for prerequisite in [tadfile, conf_file]:
		if os.path.exists(prerequisite): continue
		print(prerequisite, "not found")
		exit()

	db = connect_to_mysql(conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,'progesterone')
	qry = "select r.chromosome, r.rfrom, r.rto, r.strand from regions as r, genes as g "
	qry += "where g.name='%s' and g.region_id=r.id and r.assembly='%s' " % (gene_name,assembly)
	ret = search_db(cursor,qry)
	if not ret or (type(ret[0][0])==str and 'Error' in ret[0][0]):
		search_db(cursor,qry, verbose=True)
		exit()
	cursor.close()
	db.close()

	# there might be multiple returns, corresponding to different splices
	[chromosome, min_start, max_end, strand] = ret[0]
	for row in ret:
		[chromosome, start, end, strand] = row
		min_start = start if min_start>start else min_start
		max_end = end if max_end<end else max_end


	print ( "{} {} {}:{}-{}".format(gene_name, strand, chromosome, min_start, max_end) )

	[start, end] = get_tad (tadfile, chromosome, [min_start,max_end])
	print ("TAD containing %s region: %s:%d-%d   length %d"%(gene_name, chromosome, start, end, end-start+1))


#########################################
if __name__ == '__main__':
	main()




