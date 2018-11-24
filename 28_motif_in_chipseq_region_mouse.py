#!/usr/bin/python3


from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
import os
from utils.utils import *


#########################################
def read_tfbs_ranges(infile, tf_name):
	chipseq_regions =[]
	for line in open(infile, "r"):
		if line[0]=='%': continue
		[chrom, chromStart, chromEnd, name] = line.rstrip().split("\t")[:4]
		if name!=tf_name:continue
		chipseq_regions.append("{}_{}".format(chromStart, chromEnd))
	return chipseq_regions

#########################################
def main():

	tf = "PGR"
	species = "human"

	if species=="human":
		assembly = "hg19"
		chromosome = "4"
	else:
		assembly = "mm9"
		chromosome = "8"


	tfbs_file = "raw_data/tf_binding_sites/Hand2_{}_tfbs_{}.tsv".format(tf,assembly)
	jaspar_motifs_file  = "/storage/databases/jaspar/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt"
	if tf=="PGR":
		jaspar_motifs_file  = "/storage/databases/hocomoco/HOCOMOCOv11_core_MOUSE_mono_jaspar_format.txt"
	chipseq_regions_dir = "raw_data/chipseq_regions_%s_seqs" % species
	alignments_dir      = "raw_data/alignments_%s" % species
	for f in [tfbs_file, jaspar_motifs_file, chipseq_regions_dir]:
		if not os.path.exists(f):
			print(f,"not found")
			exit()
	if tf == "PGR":
		motif = read_pfm(jaspar_motifs_file, 'PRGR_MOUSE.H11MO.0.A')
	else:
		motif = read_pfm(jaspar_motifs_file, tf)

	# add something so that the counts are not 0
	pwm = motif.counts.normalize(pseudocounts=1)
	pssm = pwm.log_odds()




	chipseq_regions = read_tfbs_ranges(tfbs_file, tf)
	for region in chipseq_regions:
		[start, end] = [int(i) for i in region.split("_")]
		seq = read_or_download_sequence(chipseq_regions_dir, assembly, chromosome, tf, start, end)
		bpseq = Seq(seq,unambiguous_dna)
		reasonable_hits = ""
		for position, score in pssm.search(bpseq, threshold=10.0):
			revstrand=False
			if position>0:
				offset = position
				reasonable_hits += ("offset: %d score = %5.1f" % (offset, score)) + "\n"
				matched_seq = bpseq[position:position+motif.length]
			else:
				revstrand=True
				offset = len(bpseq)+position
				reasonable_hits += ("offset: %d, motif score = %5.1f (on the compl strand)" % (offset, score)) + "\n"
				matched_seq = bpseq[position:position+motif.length].reverse_complement()
			reasonable_hits += (matched_seq.upper()) + "\n"
			reasonable_hits += (motif.consensus) + "  <--- consensus"+ "\n"
			reasonable_hits += bpseq[position:position+motif.length]+ "  <--- direct strand" + "\n"
			almtfile = get_alignment_file(alignments_dir, species, assembly, chromosome,
			                                  tf, start+offset, start+offset+motif.length)
			almt =  almt_simplified(species, almtfile,pssm, revstrand)
			reasonable_hits += (almt) + "\n"

		if len(reasonable_hits)>0:
			print ("********************************")
			print ("{}, {}, {} binding site candidates".format(species, assembly, tf))
			print ("ChIPSeq region:", start, end)
			print(reasonable_hits)


#########################################
########################################
if __name__ == '__main__':
	main()
