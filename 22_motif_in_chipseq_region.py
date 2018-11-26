#!/usr/bin/python3


from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
import sys
from utils.utils import *


#########################################
def read_tfbs_ranges(infile, chr,  tf_name):
	chipseq_regions =[]
	for line in open(infile, "r"):
		if line[0]=='%': continue
		[chrom, chromStart, chromEnd, name] = line.rstrip().split("\t")[:4]
		if name!=tf_name:continue
		if chrom!="chr%s"%chr:continue
		chipseq_regions.append("{}_{}".format(chromStart, chromEnd))
	if len(chipseq_regions)==0:
		print("No chipseq region found for chr%s, transcription factor %s in %s." %(chr, tf_name, infile))
		exit()
	return chipseq_regions

#########################################
def main():

	if len(sys.argv)<6:
		print("Usage: %s <species> <gene_name> <tf_name> <tfbs_file> <'almts'|'no_almts'>" % sys.argv[0])
		exit()

	do_almts = True
	[species, gene_name, tf_name, tfbs_file, almts_choice] = sys.argv[1:6]
	if almts_choice=='no_almts': do_almts=False

	if species=="human":
		assembly = "hg19"
		chromosome = "4"
	elif species=="mouse":
		assembly = "mm9"
		chromosome = "8"
	else:
		print ("%s not enabled" % species)
		exit()


	if tf_name=="PGR":
		motifs_file  = "/storage/databases/hocomoco/HOCOMOCOv11_core_%s_mono_jaspar_format.txt" % species.upper()
	else:
		motifs_file  = "/storage/databases/jaspar/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.txt"

	chipseq_regions_dir = "raw_data/chipseq_regions_%s_seqs" % species
	alignments_dir      = "raw_data/alignments_%s" % species

	for dependency in [tfbs_file, motifs_file, chipseq_regions_dir, alignments_dir]:
		if not os.path.exists(dependency):
			print(dependency,"not found")
			exit()

	if tf_name == "PGR":
		motif = read_pfm(motifs_file, "PRGR_%s.H11MO.0.A"%species.upper())
	else:
		motif = read_pfm(motifs_file, tf_name)

	# add something so that the counts are not 0
	pwm = motif.counts.normalize(pseudocounts=1)
	pssm = pwm.log_odds()


	chipseq_regions = read_tfbs_ranges(tfbs_file, chromosome, tf_name)
	reasonable_hits_any = False
	for region in chipseq_regions:
		[start, end] = [int(i) for i in region.split("_")]
		seq = read_or_download_sequence(chipseq_regions_dir, assembly, chromosome, tf_name, start, end)
		bpseq = Seq(seq,unambiguous_dna)
		reasonable_hits = ""
		for position, score in pssm.search(bpseq, threshold=5.0):
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
			if do_almts:
				almtfile = get_alignment_file(alignments_dir, species, assembly, chromosome,
										tf_name, start+offset, start+offset+motif.length)
				almt =  almt_simplified(species, almtfile,pssm, revstrand)
				reasonable_hits += (almt) + "\n"

		if len(reasonable_hits)>0:
			reasonable_hits_any = True
			print ("********************************")
			print ("{}, {}, {} binding site candidates".format(species, assembly, tf_name))
			print ("ChIPSeq region:", start, end)
			print(reasonable_hits)

	if not reasonable_hits_any:
		print("No hits found. Try lowering the score threshold in motif search.")

#########################################
########################################
if __name__ == '__main__':
	main()
