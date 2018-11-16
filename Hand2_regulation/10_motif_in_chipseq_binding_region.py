#!/usr/bin/python3

from  Bio import motifs
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
import os
from linkto_python_modules.utils import *
import numpy as np

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
def almt_simplified(almtfile,pssm,revstrand):
	outf = open(almtfile, "r")
	assmbs = []
	seq = {}
	ranges = {}
	for line in outf:
		line = line.rstrip()
		if line[0]=='>':
			[asm,range] = line[1:].split(':')
			if not asm in seq:
				assmbs.append(asm)
				seq[asm] = ""
				ranges[asm]=[]
			ranges[asm].append(range)
		else:
			seq[asm]+= line
	outf.close()

	almt = ""
	for asm in assmbs:
		fields = asm.split(".")
		[species, chrom]  =  fields[:2]
		seq_straight = seq[asm].replace("-","")[:20]
		biopythonseq = Seq(seq_straight,unambiguous_dna)
		if revstrand:
			biopythonseq = biopythonseq.reverse_complement()
		try:
			score = pssm.calculate(biopythonseq)
			maxscore = np.amax(score)
		except:
			maxscore = -100
		almt += "%-10s %-20s %5.1f   %-15s %s\n"%(species, seq_straight, maxscore, chrom, ranges[asm])
		if (species=='rn5'):
			almt+= "-------------------------------------------\n"
	return almt

#####
def get_alignment_file(alignments_dir, species, assembly, chromosome, tf, start, end):
	almtfile = "{}/{}_{}_{}_{}_{}.txt".format(alignments_dir, tf, assembly, chromosome, start, end)
	if not (os.path.exists(almtfile) and os.path.getsize(almtfile)>0):
		get_alignment(species, assembly, chromosome, start, end, alignments_dir, almtfile)
	return almtfile

#########################################
def main():

	tf = "ESR1"
	species = "human"
	assembly = "hg19"
	chromosome = "4"

	tfbs_file = "raw_data/esr1_binding.tsv"
	jaspar_motifs_file  = "/storage/databases/jaspar/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
	chipseq_regions_dir = "raw_data/chipseq_regions"
	alignments_dir      = "raw_data/alignments"
	for f in [tfbs_file, jaspar_motifs_file, chipseq_regions_dir]:
		if not os.path.exists(f):
			print(f,"not found")
			exit()

	motif = read_pfm(jaspar_motifs_file, tf)
	# add something so that the counts are not 0
	pwm = motif.counts.normalize(pseudocounts=1)
	pssm = pwm.log_odds()

	chipseq_regions = read_tfbs_ranges(tfbs_file)
	for region in chipseq_regions:
		[start, end] = [int(i) for i in region.split("_")]
		print ("\n********************************")
		print (start, end)
		seq = read_or_download_sequence(chipseq_regions_dir, assembly, chromosome, tf, start, end)
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
			almtfile = get_alignment_file(alignments_dir, species, assembly, chromosome,
			                                  tf, start+offset, start+offset+motif.length)
			almt =  almt_simplified(almtfile,pssm, revstrand)
			print(almt)

#########################################
########################################
if __name__ == '__main__':
	main()
