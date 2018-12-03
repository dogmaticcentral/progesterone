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
from utils.utils import *
from math import log
import matplotlib.pyplot as plt

import numpy as np


def rescaled (x, offset, scale):
	return (x-offset)/scale

def atac_points(db, cursor, assembly, chromosome, tad_start, tad_end, tad_length):

	points  = []
	values  = []
	weights = []

	qry  = "select r.rfrom, r.rto, a. logfold_change, a.pval from regions as r, atac_acc_changes as a "
	qry += "where r.assembly='%s' and r.chromosome='%s' " % (assembly, chromosome)
	qry += "and r.rtype='atacseq' and r.rfrom>=%d and r.rto<=%d " % (tad_start, tad_end)
	qry += "and a.atac_region_id=r.id"
	ret = search_db(cursor,qry)
	hard_check (db, cursor, ret, qry)

	for line in ret:
		[start, end, logfold_change, pval] = line
		start = int(start)
		end = int(end)
		if end<tad_start: continue
		if start>tad_end: continue
		midpoint = (start+end)/2
		point_rescaled = rescaled(midpoint, tad_start, tad_length)
		if pval==0:
			weight = 5
		else:
			weight = -log(float(pval))
		points.append(point_rescaled)
		values.append(float(logfold_change))
		weights.append(weight)

	sorted_indices = [i for i in range(len(points))]
	sorted_indices.sort(key=lambda i: points[i])

	np_points  = np.array(points)[sorted_indices]
	np_values  = np.array(values)[sorted_indices]
	np_weights = np.array(weights)[sorted_indices]
	#for i in range(np_points.size):
	#	print("%.3f   %5.1f   %1.f" % (np_points[i], np_values[i], np_weights[i]))

	return np_points, np_values, np_weights


#########################################
def main():

	assembly  = "hg19"
	gene_name           = "Hand2"
	tad_external_exp_id = "ENCFF633ORE"
	atac_pubmed_id      = "29259032"

	conf_file = "/home/ivana/.mysql_conf"
	for prerequisite in [conf_file]:
		if os.path.exists(prerequisite): continue
		print(prerequisite, "not found")
		exit()

	######################################
	db = connect_to_mysql(conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,'progesterone')

	atac_xref_id = store_xref (cursor, 'pubmed', atac_pubmed_id)
	tad_xref_id  = get_xref_id(db,cursor,tad_external_exp_id)
	# find gene coordinates
	[chromosome, strand, gene_start, gene_end] = get_gene_coords(db,cursor,gene_name,assembly)
	# tad?
	[tad_start, tad_end] = get_tad_region(db, cursor, tad_xref_id, chromosome, gene_start, gene_end)
	tad_length = tad_end - tad_start

	# binding regions
	chipseq_regions_ESR1 = get_binding_regions_in_interval(db, cursor, assembly, chromosome,  tad_start, tad_end, 'ESR1')
	chipseq_regions_PGR  = get_binding_regions_in_interval(db, cursor, assembly, chromosome,  tad_start, tad_end, 'PGR')

	# atacseq regions (np_ are numpy arrays)
	np_atac_points, np_atac_values, np_atac_weights = atac_points(db, cursor, assembly, chromosome,  tad_start, tad_end, tad_length)

	cursor.close()
	db.close()

	######################################
	colors = []
	maxw = np.amax(np_atac_weights)
	for w in np_atac_weights:
		colors.append((0.29,0,0.5, w/maxw))

	fig, ax = plt.subplots()
	ax.axhline(0, color='black', lw=1)
	#ax.plot(np_points, np_values, 'o')
	plt.bar(np_atac_points, np_atac_values, width=0.01, color=colors)

	hand2_start = rescaled(gene_start, tad_start, tad_length)
	hand2_end   = rescaled(gene_end, tad_start, tad_length)
	ax.hlines(y=0.0, xmin=hand2_start-0.002, xmax=hand2_end+0.002, linewidth=30, color='b')

	for chipseq_reg in chipseq_regions_ESR1:
		[start, end] = [rescaled(int(p), tad_start, tad_length) for p in chipseq_reg]
		ax.hlines(y=0.0, xmin=start-0.002, xmax=end+0.002, linewidth=30, color='r')

	for chipseq_reg in chipseq_regions_PGR:
		[start, end] = [rescaled(int(p), tad_start, tad_length) for p in chipseq_reg]
		ax.hlines(y=0.0, xmin=start-0.002, xmax=end+0.002, linewidth=30, color='g')

	plt.show()

	return

#########################################
if __name__== '__main__':
	main()
