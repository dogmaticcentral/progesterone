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

	species = 'mouse'
	tf_name = 'PGR'
	verbose = True

	if not species in ['human','mouse']:
		print ("%s not enabled" % species)
		exit()

	if species=="human":
		assembly = "hg19"
		chromosome = "chr4"
		sequences_dir = "raw_data/chipseq_region_seqs_human"
		score_threshold = 10
	else:
		assembly = "mm9"
		chromosome = "chr8"
		sequences_dir = "raw_data/chipseq_region_seqs_mouse"
		score_threshold = 5

	conf_file = "/home/ivana/.mysql_conf"
	if tf_name=="PGR":
		motifs_file  = "/storage/databases/hocomoco/HOCOMOCOv11_core_%s_mono_jaspar_format.txt" % species.upper()
		pubmed_id = '29140464'
	else:
		motifs_file  = "/storage/databases/jaspar/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
		pubmed_id = '29140473'

	for dependency in [motifs_file, conf_file, sequences_dir]:
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

	#########################
	# plug in to local database
	db = connect_to_mysql(conf_file)
	cursor = db.cursor()
	search_db(cursor,"set autocommit=1")
	switch_to_db(cursor,'progesterone')

	#########################
	# find or store the reference
	xref_id = store_xref(cursor, 'pubmed', pubmed_id)

	#########################
	# read chipseq_regions
	chipseq_regions = get_binding_regions(db, cursor, assembly, chromosome, tf_name, return_binding_site_id=True)

	#########################
	# search_regions
	reasonable_hits_any = False
	for [chipseq_site_id, start, end] in chipseq_regions:
		seq = read_or_download_sequence(sequences_dir, assembly, chromosome, tf_name, start, end)
		bpseq = Seq(seq,unambiguous_dna)
		reasonable_hits = "" # for verbose output, otherwise ignored
		for position, score in pssm.search(bpseq, threshold=score_threshold):
			if position>0:
				offset = position
				matched_seq = bpseq[position:position+motif.length]
				strand = '+'
				if verbose: reasonable_hits += ("offset: %d score = %5.1f" % (offset, score)) + "\n"
			else:
				offset = len(bpseq)+position
				matched_seq = bpseq[position:position+motif.length].reverse_complement()
				strand = '-'
				if verbose:  reasonable_hits += ("offset: %d, motif score = %5.1f (on the compl strand)" % (offset, score)) + "\n"
			if verbose:
				reasonable_hits += (matched_seq.upper()) + "\n"
				reasonable_hits += (motif.consensus) + "  <--- consensus"+ "\n"
				reasonable_hits += bpseq[position:position+motif.length] + "  <--- direct strand" + "\n"

			# store region (address)
			fields  = {'species':species, 'chromosome':chromosome,'assembly':assembly, 'rtype':'motif',
							'rfrom': start+offset , 'rto':start+offset+len(matched_seq)-1, 'strand':strand, 'xref_id':xref_id}
			region_id = store_or_update(cursor, 'regions', fields, None)

			# store motif
			fixed_fields = {'region_id': region_id, 'tf_name': tf_name, 'sequence': str(matched_seq.upper()),
							'consensus': str(motif.consensus.upper()), 'score':score, 'xref_id':xref_id}
			motif_id = store_or_update(cursor, 'motifs', fixed_fields, None)

			# store mapping between the two
			fixed_fields = {'binding_site_id':chipseq_site_id, 'motif_id':motif_id}
			store_or_update(cursor, 'binding_site2motif', fixed_fields, None)

		if len(reasonable_hits)>0:
			reasonable_hits_any = True
			print ("********************************")
			print ("{}, {}, {} binding site candidates".format(species, assembly, tf_name))
			print ("ChIPSeq region:", start, end)
			print(reasonable_hits)

	if not reasonable_hits_any:
		print("No hits found. Try lowering the score threshold in motif search.")

	cursor.close()
	db.close()

#########################################
########################################
if __name__ == '__main__':
	main()
