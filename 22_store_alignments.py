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


def get_motif_regions_wo_alignment(db, cursor):
	qry  = "select m.id, r.species, r.assembly, r.chromosome, r.rfrom, r.rto, r.strand  "
	qry += "from motifs as m, regions as r "
	qry += "where m.alignment_id is NULL "
	qry += "and m.region_id=r.id"
	ret = search_db(cursor,qry)
	hard_check(db,cursor,ret,qry)
	return ret

#########################################
def main():

	conf_file  = "/home/ivana/.mysql_conf"
	scratch    = "/home/ivana/scratch"
	for dependency in [conf_file, scratch]:
		if not os.path.exists(dependency):
			print(dependency,"not found")
			exit()

	#########################
	# plug in to local database
	db = connect_to_mysql(conf_file)
	cursor = db.cursor()
	search_db(cursor,"set autocommit=1")
	switch_to_db(cursor,'progesterone')

	for address in get_motif_regions_wo_alignment(db, cursor):
		print(address)
		[motif_id, species, assembly, chrom, region_from, region_to, strand] = address
		almtfile = "{}/{}.afa".format(scratch, os.getpid())
		get_alignment(species, assembly, chrom, region_from, region_to, scratch, almtfile)
		print(almtfile)
		names_ordered, almt = almt_simplified(almtfile, species)
		for nm,sequence in almt.items():
			print(">"+nm)
			print(sequence)

		# store each sequence from this alignment as motif
		motif_ids = []
		sequences = []
		for name in names_ordered:
			[asm, chromrange] = name.split(".")
			[chrom, range] = chromrange.split(":")
			seq_straight = almt[name].replace("-", "").upper()
			biopythonseq = Seq(seq_straight, unambiguous_dna)
			if strand=='-': biopythonseq = biopythonseq.reverse_complement()
			if asm==assembly:
				mi = motif_id
			else:
				#range_id = store_range(cursor, asm, chrom, range) # ranges table
				#mi = store_motif(cursor, range_id, str(biopythonseq)) # motifs table
				mi = -1
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
