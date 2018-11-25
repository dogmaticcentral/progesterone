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


from utils.utils import *


def read_pfm(jaspar_motifs_file, tf_name):
	motif = None
	with open(jaspar_motifs_file) as handle:
		for m in motifs.parse(handle,"jaspar"):
			if m.name==tf_name:
				motif = m
				break
	return motif

#########################################
def read_or_download_sequence(chipseq_regions_dir, assembly, chromosome, tf, start, end):
	seqfile = "{}/{}_{}_{}_{}_{}.txt".format(chipseq_regions_dir,tf, assembly, chromosome, start, end)
	if (os.path.exists(seqfile) and os.path.getsize(seqfile)>0):
		outf = open(seqfile, "r")
		seq = outf.read()
		outf.close()
	else:
		seq = ucsc_fragment_sequence(assembly, chromosome, start, end)
		outf = open(seqfile, "w")
		outf.write(seq.replace("\n",""))
		outf.close()

	return seq


#########################################
# noinspection PyUnreachableCode
def main():

	tf = "ESR1"

	if True:
		assembly = "hg19"
		chromosome = "4"
		start = 174447651- 1000000
		end   = 174447651
	else:
		#assembly = "mm10"
		#chromosome = "8"
		#Hand2 at chr8:57320983-57324517
		#start = 57324517
		#end = 57324517 + 1000000
		assembly = "rn4"
		# Hand2 at chr16:36324327-36326037
		chromosome = "16"
		start = 36326037
		end = 36326037 + 1000000

	jaspar_motifs_file  = "/storage/databases/jaspar/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
	for f in [ jaspar_motifs_file]:
		if not os.path.exists(f):
			print(f,"not found")
			exit()

	motif = read_pfm(jaspar_motifs_file, tf)
	# add something so that the counts are not 0
	pwm = motif.counts.normalize(pseudocounts=1)
	pssm = pwm.log_odds()

	# seq from UCSC
	print("intial range:", start, end)
	seq = ucsc_fragment_sequence(assembly,chromosome,start,end)
	bpseq = Seq(seq,unambiguous_dna)

	for position, score in pssm.search(bpseq, threshold=10.0):
		if position>0:
			offset = position
			print("range:", start, end)
			print("offset %d: score = %5.1f" % (offset, score))
			matched_seq = bpseq[position:position+motif.length]
		else:
			offset = len(bpseq)+position
			print("offset %d: score = %5.1f (on the compl strand)" % (offset, score))
			matched_seq = bpseq[position:position+motif.length].reverse_complement()
		print(motif.consensus)
		print(matched_seq.upper())
		print(bpseq[position:position+motif.length],"  <--- direct strand")
		print()

#########################################
########################################
if __name__ == '__main__':
	main()
