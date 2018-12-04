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
	report_vals = []
	selfint_value = find_selfint(db, cursor, assembly, chromosome, gene_name)
	if not selfint_value:
		return
	ints = get_interacting_regions(db, cursor, assembly, chromosome, gene_name, tfbs_start, tfbs_end)
	for [rfrom, rto, interaction] in ints:
		#print("\tinteraction:", rfrom, rto, interaction, "%.1f%%"%(interaction/selfint_value*100))
		address = "{}:{}:{}-{}".format(assembly, chromosome, rfrom, rto)
		report_vals.append([address, interaction, "%.1f%%"%(interaction/selfint_value*100)])
	return report_vals


#########################################
def find_dist_to_gene (strand, gene_start, gene_end, rfrom, rto):
	if rto<gene_start:
		dist = (gene_start-rto)/1000
		direction = "upstream" if strand=="+" else "downstream"
	elif gene_end<rfrom:
		dist = (rfrom-gene_end)/1000
		direction = "downstream" if strand=="+" else "upstream"
	else:
		dist = 0
		direction = "within gene_region"

	return "%1.fKbp %s" %(dist, direction)


#########################################
def get_alignment(db, cursor, alignment_id, gene_name, consensus_score):
	qry = "select motif_ids, alignment from alignments where id=%d" % alignment_id
	ret = search_db(cursor,qry)
	hard_check(db, cursor, ret, qry)
	motif_ids = ret[0][0].split(",")
	seqs = ret[0][1].split(",")
	seq = dict(zip(motif_ids,seqs))
	# get info for all motifs
	qry  = "select m.id, m.consensus, m.score, r.assembly, r.chromosome, r.rfrom, r.rto, r.strand "
	qry += "from regions as r, motifs as m "
	qry += "where m.alignment_id =%d  and m.region_id=r.id" % alignment_id
	ret = search_db(cursor,qry)
	hard_check(db, cursor, ret, qry)
	label = {}
	score = {}
	label['0'] = "consensus"
	seq['0'] = None
	score['0'] = "%5.1f"%consensus_score
	for [motif_id, consensus, scr,  assembly, rchrom, rfrom, rto, rstrand]  in ret:
		if not seq['0']: seq['0'] = consensus
		# species common name
		species_common = assembly2species_common(cursor,assembly)
		# gene position in this assembly
		[gene_chr, gene_strand, gene_from, gene_to] = get_gene_coords (db, cursor, gene_name, assembly)
		# distance to gene
		if rchrom!=gene_chr:
			distance_to_gene = "(? contig)"
		else:
			distance_to_gene = find_dist_to_gene(gene_strand, gene_from, gene_to, rfrom, rto)
		address = "{}:{}:{}-{}:{}".format(assembly,rchrom,rfrom,rto,rstrand)
		label[str(motif_id)] = ", ".join([species_common, address,distance_to_gene])
		score[str(motif_id)] = "%5.1f"%scr
	# motif ids are in the order in which we want to see the alignment
	return [['0']+motif_ids, seq, label, score]


#########################################
def alignment_report(motif_ids, seq, label, score):
	for mi in motif_ids:
		if '-' in seq[mi]: return ""
	ret = "----------------\n"
	for mi in motif_ids:
		ret += "%-62s %20s %5s (<--motif score)\n"%(label[mi], seq[mi], score[mi])
	return ret


#########################################
def report (db, cursor, assembly,  species, gene_name, tf_name, tad_external_exp_id):

	if tf_name=="PGR":
		motifs_file  = "/storage/databases/hocomoco/HOCOMOCOv11_core_%s_mono_jaspar_format.txt" % species.upper()
		motif = read_pfm(motifs_file, "PRGR_%s.H11MO.0.A"%species.upper())
	else:
		motifs_file  = "/storage/databases/jaspar/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
		motif = read_pfm(motifs_file, tf_name)
	# add something so that the counts are not 0
	pwm = motif.counts.normalize(pseudocounts=1)
	pssm = pwm.log_odds()
	consensus = motif.consensus
	cons_score = pssm.calculate(consensus)

	# find gene coordinates
	[chromosome, gene_strand, gene_start, gene_end] = get_gene_coords(db,cursor,gene_name,assembly)
	# tad?
	if species=='human':
		tad_xref_id  = get_xref_id(db,cursor, tad_external_exp_id)
		[tad_start, tad_end] = get_tad_region(db, cursor, tad_xref_id, chromosome, gene_start, gene_end)
	else:
		[tad_start, tad_end] = [gene_start - 1.e6, gene_end+ 1.e6]

	# get all tfbs' inside tad
	tfbs = get_binding_regions_in_interval(db, cursor, assembly,
											chromosome, tad_start, tad_end, tf_name, return_binding_site_id=True)
	# for each tfbs
	for [binding_site_id, tfbs_start, tfbs_end] in tfbs:

		tfbs_address   = "{}:{}:{}-{}".format(assembly,chromosome,tfbs_start, tfbs_end)
		tfbs_distance  = find_dist_to_gene (gene_strand, gene_start, gene_end, tfbs_start, tfbs_end)
		qry = "select xref_id from binding_sites where id=%d" % binding_site_id
		ret = search_db(cursor, qry)
		hard_check(db, cursor, ret, qry)
		xref_id= ret[0][0]
		tfbs_references = get_references(db, cursor, xref_id)
		#   find interacting region(s) and interaction value(s)
		interaction_rept = int_report(db, cursor, assembly, chromosome, gene_name, tfbs_start, tfbs_end)
		# find nearby accessibility change regions
		# for atac_region  in get_atac_regions(db, cursor, assembly, chromosome, tfbs_start-10000, tfbs_end+10000):
		# 	print ("\tatac region:", atac_region)
		#   find all motifs
		motif_ids = get_motifs_in_binding_site(db, cursor, binding_site_id)
		#   for each motif
		almt = ""
		for motif_id in motif_ids:
			# find motif score, alignment
			qry = "select * from motifs where id=%d" % motif_id
			ret = search_db(cursor, qry)
			hard_check(db,cursor, ret, qry)
			[mid, region_id, tf_name, sequence, cons, score, xref_id, alignment_id] = ret[0]
			if cons != consensus:
				print("something's wrong: consensus mismatch", cons, consensus)
				exit()
			# reconstruct alignment from alignment_id
			[motif_ids, seqs, labels, score] = get_alignment(db, cursor, alignment_id, gene_name, cons_score)
			almt += alignment_report(motif_ids, seqs, labels, score)
		# this could be formatted into a proper table or some such
		print("-----------------------------")
		print("ChipSeq binding site: {} {}".format(tfbs_address, tfbs_distance))
		print("Sources(s): {}".format(tfbs_references.replace(",ENCODE:",",")))
		if interaction_rept:
			for [interval_addres, ivalue, pct] in interaction_rept:
				print("Hic interaction interval: {}, {} of self-interaction".format(interval_addres, pct))
		if len(almt)>0: print(almt)
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
			print("====================================")
			print ("\n\n{} {}".format(species, tf_name))
			report (db, cursor, assembly[species], species, gene_name, tf_name, tad_external_exp_id)

	cursor.close()
	db.close()


#########################################
########################################
if __name__ == '__main__':
	main()
