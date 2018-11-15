#!/usr/bin/python3

from  Bio import motifs
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
import os
from linkto_python_modules.utils import *


#########################################
def read_tfbs_ranges(infile):
	chipseq_regions =[]
	for line in open(infile, "r"):
		if line[0]=='%': continue
		[name, chipseq_region, encode_tfbs_score, hic_region, hic_int_score] = line.rstrip().split("\t")
		chipseq_regions.append(chipseq_region)
	return chipseq_regions

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
def main():

	tf = "ESR1"
	assembly = "hg19"
	chromosome = "4"

	tfbs_file = "raw_data/esr1_binding.tsv"
	jaspar_motifs_file = "/storage/databases/jaspar/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
	chipseq_regions_dir = "raw_data/chipseq_regions"
	for f in [tfbs_file, jaspar_motifs_file, chipseq_regions_dir]:
		if not os.path.exists(f):
			print(f,"not found")
			exit()

	motif = read_pfm(jaspar_motifs_file, 'ESR1')
	# add something so that the counts are not 0
	pwm = motif.counts.normalize(pseudocounts=1)
	pssm = pwm.log_odds()

	chipseq_regions = read_tfbs_ranges(tfbs_file)
	for region in chipseq_regions:
		[start, end] = [int(i) for i in region.split("_")]
		print ()
		print (start, end)
		seq = read_or_download_sequence(chipseq_regions_dir, assembly, chromosome, tf, start, end)
		bpseq = Seq(seq,unambiguous_dna)
		for position, score in pssm.search(bpseq, threshold=3.0):
			print("Position %d: score = %5.3f" % (position, score))
			if position>0:
				matched_seq = bpseq[position:position+motif.length]
			else:
				matched_seq = bpseq[position:position+motif.length].reverse_complement()
			print(motif.consensus)
			print(matched_seq.upper())


#########################################
########################################
if __name__ == '__main__':
	main()
