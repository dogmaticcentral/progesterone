#!/usr/bin/python3

# source https://www.encodeproject.org/experiments/ENCSR551IPY/
# hdf files can be inspected with h5dump --contents <filename>
# even easier: h5ls <filename>  (or h5ls -vlr <filename> but this might be too verbose)
import math

import h5py
import os
from utils.utils import *

def range2tfb(infile):
	r2name  = {}
	chrom = None
	for line in open(infile, "r"):
		[chrom_readin, chromStart, chromEnd, name] = line.rstrip().split("\t")[:4]
		if chrom==None:
			chrom=chrom_readin
		elif (chrom != chrom_readin):
			print ("Unexpected: chrom not the same for all entry lines in %s", infile)
			exit()
		range = "{}_{}".format(chromStart,chromEnd)
		if not range in r2name:
			r2name[range]=[]
		r2name[range].append(name)

	return r2name


def inside(container_from, container_to, start, end):
	return container_from<=start<=container_to or container_from<=end<=container_to


def find_tfbs_in(bfrom, bto, r2name):
	tfbs = []
	for range, names in r2name.items():
		[start, end]= [int(i) for i in range.split("_")]
		if inside(bfrom, bto, start, end):
			tfbs.append(names +  [range])
	return tfbs


#########################################
def main():
	# from 02_emve_tads.py
	gene_name = "Hand2"
	tfbs_dir  = "raw_data/tf_binding_sites"
	tfbs_file = "raw_data/tf_binding_sites/Hand2_ESR1_tfbs_hg19.tsv"

	ucsc_gene_regions_dir = "/storage/databases/ucsc/gene_ranges/human/hg19"
	tad_file = "/storage/databases/tads/encode/ENCFF633ORE.bed"
	contacts_file= "/storage/databases/encode/ENCSR551IPY/ENCFF331ABN.h5"
	hic_dir = "raw_data/hic_interactions"

	for f in [ucsc_gene_regions_dir,tad_file, tfbs_dir, tfbs_file, contacts_file, hic_dir]:
		if not os.path.exists(f):
			print(f,"not found")
			exit()

	chromosome, strand, gene_range = ucsc_gene_coords(gene_name, ucsc_gene_regions_dir)
	if not chromosome or not gene_range:
		print("coordinates not found for", gene_name)
		exit()

	promoter_range = []
	if strand=='+':
		promoter_range = [gene_range[0]-1000, gene_range[0]]
	elif strand=='-':
		promoter_range = [gene_range[1], gene_range[1]+1000]
	else:
		print ("strand %s not recognized"%strand)
		exit()

	print("%s  chromsome: %s   promoter range:"%(gene_name, chromosome), promoter_range)


	[tad_start, tad_end] = get_tad (tad_file, chromosome, gene_range)

	infile= h5py.File(contacts_file,'r')
	chrom_index = [i for i in range(24) if infile['chrs'][i]==bytes(chromosome,'ascii')][0]
	print("%s  chromsome index in the hdf5 file: %d"% (chromosome, chrom_index))
	chrom_bin_from, chrom_bin_to = infile['chr_bin_range'][chrom_index]
	bin_positions = infile['bin_positions']
	pstart = promoter_range[0]
	pend   = promoter_range[1]
	#promoter_bins = [b for b in range(bin_from,bin_to+1) if
	#				(bin_positions[b][1]<=pstart<=bin_positions[b][2]) or
	#				(bin_positions[b][1]<=pend<=bin_positions[b][2])]
	# this is 2.5 times faster
	promoter_bins = []
	for b in range(chrom_bin_from,chrom_bin_to+1):
		bfrom, bto =bin_positions[b][1:3]
		if bfrom<=pstart<=bto or  bfrom<=pend<=bto: promoter_bins.append(b)

	if len(promoter_bins)>1:
		print("bins containining promoter:", promoter_bins)
		for pb in promoter_bins:
			print (pb, infile['bin_positions'][pb])
		print ("promoter straddles two intervals; please generalize me accordingly")
		exit()

	pb = promoter_bins[0]
	print("promoter bin: ", pb)
	pb_interactions = infile['interactions'][pb][:]
	int_strength = {}
	for b in range(chrom_bin_from, chrom_bin_to+1):
		if not math.isnan(pb_interactions[b]) and pb_interactions[b]>0.0:
			int_strength[b]= pb_interactions[b]

	bins_sorted_by_strength = [b for b in sorted(int_strength, key=int_strength.get, reverse=True)]

	# regions we are interested in (expected usage: from the previous script, 03_tf_binding_sites_from_UCSC.py)
	r2name = range2tfb(tfbs_file)

	print()
	print ("number of bins in chromosome %s"%chromosome, chrom_bin_to-chrom_bin_from+1, "interacting bins:",  len(bins_sorted_by_strength))

	hic_file = tfbs_file.replace(tfbs_dir,hic_dir)
	outf = open(hic_file,"w")
	outf.write("\t".join(["% name", "chipseq_region", "hic_region", "hic_int_score"])+ "\n")
	rank = 0
	for b in bins_sorted_by_strength:
		bfrom, bto =bin_positions[b][1:3]
		if tad_start<=bfrom<=tad_end  or  tad_start<=bto<=tad_end:
			rank += 1
			tfbs_inside = find_tfbs_in(bfrom, bto,r2name)

			for t in tfbs_inside:
				outf.write( "\t".join( t+["{}_{}".format(bfrom, bto), "%5.0f"%(int_strength[b])]) + "\n")


	outf.close()

	return True


#########################################
########################################
if __name__ == '__main__':
	main()
