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
def store_region(cursor, species, assembly, chrom, rfrom, rto, strand, xref_id):
	fields  = {'species':species, 'chromosome':chrom, 'assembly':assembly, 'rtype':'motif',
					'rfrom':rfrom, 'rto':rto, 'strand':strand,   'xref_id':xref_id}
	region_id = store_or_update(cursor, 'regions', fields, None)
	return region_id


#########################################
def store_inferred_motif(cursor, region_id, tf_name, biopythonseq, consensus, pssm, xref_id):

	try:
		score = pssm.calculate(biopythonseq)
		maxscore = np.amax(score)
	except:
		maxscore = -100

	fixed_fields = {'region_id': region_id, 'tf_name': tf_name, 'sequence': str(biopythonseq),'xref_id':xref_id}
	update_fields = {'consensus': str(consensus), 'score':maxscore}
	motif_id = store_or_update(cursor, 'motifs', fixed_fields, update_fields)
	return motif_id


#########################################
def simplify_alignment(cursor, labels, seqs, species):
	if (species=='human'):
		relatives = ['human', 'crab-eating macaque', 'green monkey'] # macFas & chlSab
	else:
		relatives = ['mouse', 'rat']

	seq_relatives = {}
	labels_relatives = []
	for label in labels:
		assembly = label.split("_")[0]
		tgt_species = assembly2species_common(cursor,assembly)
		if not tgt_species: continue
		if not tgt_species in relatives: continue
		labels_relatives.append(label)
		seq_relatives[label] = seqs[label]

	remove_all_gaps(seq_relatives)
	return labels_relatives, seq_relatives

#########################################
def main():

	species = 'mouse'
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
		[motif_id, qry_assembly, chrom, region_from, region_to, qry_strand] = address
		labels, seqs = get_alignment(species, qry_assembly, chrom, region_from, region_to, scratch)
		labels, seqs = simplify_alignment(cursor, labels, seqs, species)

		# store each sequence from this alignment as motif
		motif_ids = []
		sequences = []
		for label in labels:
			[tgt_assembly, chrom, rfrom, rto, tgt_strand]  = label.split("_")
			rfrom = int(rfrom)
			rto   = int(rto)
			seq_straight = seqs[label].replace("-", "").upper()
			#biopythonseq = Seq(seq_straight, unambiguous_dna)
			# we search with region on "+" strand, even if the motif is on "-"
			# thus we need to take reverse compl whenever the strands are different
			#if tgt_strand!=qry_strand: biopythonseq = biopythonseq.reverse_complement()
			if tgt_assembly==qry_assembly:
				mi = motif_id
			else:
				tgt_species = assembly2species_common(cursor,tgt_assembly)
				print(tgt_assembly, tgt_species, chrom, rfrom, rto, xref_id)
				region_id = store_region(cursor, tgt_species, tgt_assembly, chrom,
										rfrom, rto, tgt_strand, xref_id) # regions table
				# we do not need to take reverse complement, because mafsInRegion already
				# works with the complement if needed (and reports tha it was found on "-" strand)
				mi = store_inferred_motif(cursor, region_id, tf_name, seq_straight,
										motif.consensus, pssm, xref_id) # motifs table
			motif_ids.append(str(mi))
			sequences.append(seqs[label]) # these might still have  gaps

		# take motif_ids and seqs (possibly with gaps) and store them in alignments table
		print(",".join(motif_ids), ",".join(sequences))

		# store alignment
		
		exit()

	cursor.close()
	db.close()


#########################################
########################################
if __name__ == '__main__':
	main()
