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


def get_atac_regions(db, cursor, assembly, chromosome, tfbs_start, tfbs_end):
	qry  = "select r.rfrom, r.rto, a. logfold_change, a.pval from regions as r, atac_acc_changes as a "
	qry += "where r.assembly='%s' and r.chromosome='%s' " % (assembly, chromosome)
	qry += "and r.rtype='atacseq' "
	qry += "and not (r.rto<%d or  r.rfrom>%d) " % (tfbs_start, tfbs_end)
	qry += "and a.atac_region_id=r.id"
	ret = search_db(cursor,qry)
	permissive_check (db,cursor, ret, qry)
	return ret if ret else []


#########################################
def get_interacting_regions(db, cursor, assembly, chromosome, gene_name, tfbs_start, tfbs_end):
	qry  = "select r.rfrom, r.rto, h.interaction from regions as r, hic_interactions as h "
	qry += "where r.assembly='%s' and r.chromosome='%s' and r.rtype='interacting' " % (assembly, chromosome)
	qry += "and h.gene_name='%s' and  h.interacting_hic_region_id=r.id  " % (gene_name)
	qry += "and not (r.rto<%d or  r.rfrom>%d)" % (tfbs_start, tfbs_end)
	ret = search_db(cursor,qry)
	permissive_check (db,cursor, ret, qry)
	return ret if ret else []


#########################################
def find_selfint(db, cursor, assembly, chromosome, gene_name):
	qry  = "select h.interaction from regions as r, hic_interactions as h "
	qry += "where r.assembly='%s' and r.chromosome='%s' and r.rtype='interacting' " % (assembly, chromosome)
	qry += "and h.gene_name='%s' and  h.interacting_hic_region_id=r.id  " % (gene_name)
	qry += "and  h.interacting_hic_region_id = h.gene_hic_region_id"
	ret = search_db(cursor,qry)
	if not ret: return None
	permissive_check (db,cursor, ret, qry)
	if len(ret)>1:
		print(" ? multiple self int for ", assembly, chromosome, gene_name)
		cursor.close()
		db.close()
		exit()
	if ret[0][0]==0:
		print(" ? zero self int for ", assembly, chromosome, gene_name)
		cursor.close()
		db.close()
		exit()
	return ret[0][0]


#########################################
def int_report(db, cursor, assembly, chromosome, gene_name, tfbs_start, tfbs_end):
	selfint_value = find_selfint(db, cursor, assembly, chromosome, gene_name)
	if not selfint_value:
		return
	ints = get_interacting_regions(db, cursor, assembly, chromosome, gene_name, tfbs_start, tfbs_end)
	for [rfrom, rto, interaction] in ints:
		print("\t interaction:", rfrom, rto, interaction, "%.1f%%"%(interaction/selfint_value*100))


def find_dist_to_gene (strand, gene_start, gene_end,rfrom,rto):
	if rto<gene_start:
		dist = (gene_start-rto)/1000
		direction = "upstream" if strand=="+" else "downstream"
	elif gene_end<rfrom:
		dist = (rfrom-gene_end)/1000
		direction = "downstream" if strand=="+" else "upstream"
	else:
		dist = 0
		direction = "within gene_region"

	return "%1.f Kbp, %s" %(dist, direction)


#########################################
def report (db, cursor, assembly,  species, gene_name, tf_name, tad_external_exp_id):
	# find gene coordinates
	[chromosome, gene_strand, gene_start, gene_end] = get_gene_coords(db,cursor,gene_name,assembly)
	# tad?
	if species=='human':
		tad_xref_id  = get_xref_id(db,cursor, tad_external_exp_id)
		[tad_start, tad_end] = get_tad_region(db, cursor, tad_xref_id, chromosome, gene_start, gene_end)
	else:
		[tad_start, tad_end] = [gene_start - 1.e6, gene_end+ 1.e6]

	# get all tfbs' inside tad
	tfbs = get_binding_regions_in_interval(db, cursor, assembly, chromosome, tad_start, tad_end, tf_name, return_binding_site_id=True)
	# for each tfbs
	for [binding_site_id, tfbs_start, tfbs_end] in tfbs:
		print(tfbs_start, tfbs_end)
		#   find interacting region(s) and interaction value(s)
		int_report(db, cursor, assembly, chromosome, gene_name, tfbs_start, tfbs_end)
		# find nearby accessibility change regions
		for atac_region  in get_atac_regions(db, cursor, assembly, chromosome, tfbs_start-10000, tfbs_end+10000):
			print ("\tatac region:", atac_region)
		#   find all motifs
		motif_ids =  get_motifs_in_binding_site(db, cursor, binding_site_id)
		#   for each motif
		for motif_id in motif_ids:
			print("\t-------------------")
			print("\t motif_id:",motif_id)
			# find motif score, alignment
			qry = "select * from motifs where id=%d" % motif_id
			ret = search_db(cursor, qry)
			hard_check(db,cursor, ret, qry)
			[mid, region_id, tf_name, sequence, consensus, score, xref_id, alignment_id] = ret[0]
			# reconstruct alignment from alignment_id
			# find distance to gene
			[chromosome, rfrom, rto, strand] = get_region_coords (db, cursor, region_id)
			dist_to_gene = find_dist_to_gene(gene_strand, gene_start, gene_end,rfrom,rto)
			print ("\t{} ({})  score: {}  dist to gene: {}".format(sequence, consensus,score, dist_to_gene))
		#exit()
	return

#########################################
def main():

	assembly   = {'human':'hg19', 'mouse':'mm9'}
	gene_name  = "Hand2"
	tad_external_exp_id = "ENCFF633ORE"
	conf_file = "/home/ivana/.mysql_conf"
	for prerequisite in [conf_file]:
		if os.path.exists(prerequisite): continue
		print(prerequisite, "not found")
		exit()

	db = connect_to_mysql(conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,'progesterone')

	for species in ['human','mouse']:
		for tf_name in ['ESR1', 'PGR']:
			print ("\n\n", species,tf_name)
			report (db, cursor, assembly[species], species, gene_name, tf_name, tad_external_exp_id)

	cursor.close()
	db.close()


#########################################
########################################
if __name__ == '__main__':
	main()
