#!/usr/bin/python3

from  Bio import motifs
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
import os
from linkto_python_modules.utils import *
import numpy as np


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

	if False:
		assembly = "hg19"
		chromosome = "4"
		start = 174447651- 1000000
		end   = 174447651
	else:
		#assembly = "mm10"
		#chromosome = "8"
		#Hand2 at chr8:57320983-57324517
		assembly   = "calJac3"
		chromosome = "3"
		#HAND2 at chr3:15951910-15954607
		start = 15954607
		end   = 15954607 + 1000000



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
	seq = ucsc_fragment_sequence(assembly,chromosome,start,end)
	bpseq = Seq(seq,unambiguous_dna)

	for position, score in pssm.search(bpseq, threshold=10.0):
		revstrand=False
		if position>0:
			offset = position
			print("offset %d: score = %5.1f" % (offset, score))
			matched_seq = bpseq[position:position+motif.length]
		else:
			revstrand=True
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
