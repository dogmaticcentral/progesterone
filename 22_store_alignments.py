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
def get_motif_regions_wo_alignment(db, cursor, tf_name, species):
	qry  = "select m.id, r.assembly, r.chromosome, r.rfrom, r.rto, r.strand  "
	qry += "from motifs as m, regions as r "
	qry += "where m.alignment_id is NULL and m.tf_name='%s' " % tf_name
	qry += "and r.species='%s' " % species
	qry += "and m.region_id=r.id"
	ret = search_db(cursor,qry)
	hard_check(db,cursor,ret,qry)
	return ret


#########################################
def store_region(cursor, species, assembly, chrom, region, xref_id):
	[start,end] = region.split("-")
	fields  = {'species':species, 'chromosome':chrom, 'assembly':assembly, 'rtype':'motif',
					'rfrom':start, 'rto':end, 'xref_id':xref_id}
	region_id = store_or_update(cursor, 'regions', fields, None)
	return region_id


#########################################
def store_inferred_motif(cursor, region_id, tf_name, biopythonseq, consensus, pssm, xref_id):

	try:
		score = pssm.calculate(biopythonseq)
		maxscore = np.amax(score)
	except:
		maxscore = -100

	fixed_fields = {'region_id': region_id, 'tf_name': tf_name, 'sequence': str(biopythonseq),
							'consensus': consensus, 'score':maxscore, 'xref_id':xref_id}
	motif_id = store_or_update(cursor, 'motifs', fixed_fields, None)
	return motif_id

#########################################
def main():

	species = 'human'
	tf_name = 'PGR'


	conf_file  = "/home/ivana/.mysql_conf"
	scratch    = "/home/ivana/scratch"
	if tf_name=="PGR":
		motifs_file  = "/storage/databases/hocomoco/HOCOMOCOv11_core_%s_mono_jaspar_format.txt" % species.upper()
		pwm_pubmed_id = '29140464'
	else:
		motifs_file  = "/storage/databases/jaspar/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
		pwm_pubmed_id = '29140473'
	for dependency in [conf_file, scratch, motifs_file]:
		if not os.path.exists(dependency):
			print(dependency,"not found")
			exit()

	#########################
	# read/normalize PWM
	if tf_name == "PGR":
		motif = read_pfm(motifs_file, "PRGR_%s.H11MO.0.A"%species.upper())
	else:
		motif = read_pfm(motifs_file, tf_name)
	# add something so that the counts are not 0
	pwm = motif.counts.normalize(pseudocounts=1)
	pssm = pwm.log_odds()

	# multiz reference for the maf alignments downloaded from UCSC
	maf_pubmed_id = '15060014'

	#########################
	# plug in to local database
	db = connect_to_mysql(conf_file)
	cursor = db.cursor()
	search_db(cursor,"set autocommit=1")
	switch_to_db(cursor,'progesterone')

	#########################
	# store reference info
	pwm_xref_id = store_xref(cursor, 'pubmed', pwm_pubmed_id)
	maf_xref_id = store_xref(cursor, 'pubmed', maf_pubmed_id)
	# we will use combined reference for the new motifs we store
	xref_id = store_xref(cursor, 'this_db', ",".join([str(pwm_xref_id), str(maf_xref_id)]))

	for address in get_motif_regions_wo_alignment(db, cursor, tf_name, species):
		print(address)
		[motif_id, qry_assembly, chrom, region_from, region_to, strand] = address
		almtfile = "{}/{}.afa".format(scratch, os.getpid())
		get_alignment(species, qry_assembly, chrom, region_from, region_to, scratch, almtfile)
		print(almtfile)
		names_ordered, almt = almt_simplified(almtfile, species)
		for nm,sequence in almt.items():
			print(">"+nm)
			print(sequence)

		# store each sequence from this alignment as motif
		motif_ids = []
		sequences = []
		for name in names_ordered:
			[tgt_assembly, chromregion] = name.split(".")
			[chrom, region] = chromregion.split(":")
			seq_straight = almt[name].replace("-", "").upper()
			biopythonseq = Seq(seq_straight, unambiguous_dna)
			if strand=='-': biopythonseq = biopythonseq.reverse_complement()
			if tgt_assembly==qry_assembly:
				mi = motif_id
			else:
				print(tgt_assembly, chrom, region, xref_id)
				exit()
				region_id = store_region(cursor, species, tgt_assembly, chrom, region, xref_id) # regions table
				mi = store_inferred_motif(cursor, region_id, tf_name, biopythonseq,
											motif.consensus, pssm, xref_id) # motifs table
			motif_ids.append(str(mi))
			sequences.append(almt[name])

		# take motif_ids and seqs (possibly with gaps) and store them in alignments table
		print(",".join(motif_ids), ",".join(sequences))
		exit()

	cursor.close()
	db.close()


#########################################
########################################
if __name__ == '__main__':
	main()
