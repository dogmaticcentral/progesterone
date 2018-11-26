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

# source https://www.encodeproject.org/experiments/ENCSR551IPY/
# hdf files can be inspected with h5dump --contents <filename>
# even easier: h5ls <filename>  (or h5ls -vlr <filename> but this might be too verbose)
import math

import h5py
import os
from utils.utils import *

def read_ranges(infile, tf_name):
	ranges = []
	chrom = None
	for line in open(infile, "r"):
		[chrom_readin, chromStart, chromEnd, name] = line.rstrip().split("\t")[:4]
		if name != tf_name: continue
		if chrom==None:
			chrom=chrom_readin
		elif (chrom != chrom_readin):
			print ("Unexpected: chrom not the same for all entry lines in %s", infile)
			exit()
		range = "{}_{}".format(chromStart,chromEnd)
		ranges.append(range)
	return ranges


def inside(container_from, container_to, start, end):
	return container_from<=start<=container_to or container_from<=end<=container_to


def find_tfbs_in(bfrom, bto, ranges):
	tfbs = []
	for range in ranges:
		[start, end]= [int(i) for i in range.split("_")]
		if inside(bfrom, bto, start, end):
			tfbs.append(range)
	return tfbs


#########################################
# noinspection PyUnreachableCode
def main():

	assembly  = "hg19"
	gene_name = "Hand2"
	# some groseness ...
	if False:
		tf_name   = "ESR1"
		tfbs_dir  = "raw_data/tf_binding_sites_ucsc"
		tfbs_file = "{}/{}_tfbs_{}.tsv".format(tfbs_dir, gene_name, assembly)
	else:
		tf_name   = "PGR"
		tfbs_dir  = "raw_data/tf_binding_sites_geo"
		tfbs_file = "{}/{}_{}_tfbs_{}.tsv".format(tfbs_dir, gene_name, tf_name, assembly)

	ucsc_gene_regions_dir = "/storage/databases/ucsc/gene_ranges/human/hg19"
	tad_file = "/storage/databases/encode/ENCSR551IPY/ENCFF633ORE.bed"
	contacts_file = "/storage/databases/encode/ENCSR551IPY/ENCFF331ABN.h5"
	hic_dir = "raw_data/hic_interactions"

	for f in [ucsc_gene_regions_dir,tad_file, tfbs_dir, tfbs_file, contacts_file, hic_dir]:
		if not os.path.exists(f):
			print(f,"not found")
			exit()

	chromosome, strand, gene_range = ucsc_gene_coords(gene_name, ucsc_gene_regions_dir)
	if not chromosome or not gene_range:
		print("coordinates not found for", gene_name)
		exit()

	###############################
	# for promoter range take 1Kbp from the tx start region
	# it doesn't really matter, because the resolution here
	# is determined by HiC resolution, which is 40 Kbp
	promoter_range = []
	if strand=='+':
		promoter_range = [gene_range[0]-1000, gene_range[0]]
	elif strand=='-':
		promoter_range = [gene_range[1], gene_range[1]+1000]
	else:
		print ("strand %s not recognized"%strand)
		exit()

	print("%s  chromosome: %s   promoter range:" % (gene_name, chromosome), promoter_range)

	###############################
	# we are only interested in the region on the chromosome  correspoding to 'our' TAD
	[tad_start, tad_end] = get_tad (tad_file, chromosome, gene_range)

	###############################
	# HiC bins corresponding to our chromosome
	infile= h5py.File(contacts_file,'r')
	chrom_index = [i for i in range(24) if infile['chrs'][i]==bytes(chromosome,'ascii')][0]
	print("%s  chromsome index in the hdf5 file: %d"% (chromosome, chrom_index))
	chrom_bin_from, chrom_bin_to = infile['chr_bin_range'][chrom_index]

	###############################
	# all bin positions handle
	bin_positions = infile['bin_positions']

	###############################
	# bin(s) that contain the gene we are interested in
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

	###############################
	# interactions between the promoter-containing bin and everybody else
	pb_interactions = infile['interactions'][pb][:]
	int_strength = {}
	for b in range(chrom_bin_from, chrom_bin_to+1):
		if not math.isnan(pb_interactions[b]) and pb_interactions[b]>0.0:
			int_strength[b] = pb_interactions[b]
	print("number of bins in chromosome %s"%chromosome, chrom_bin_to-chrom_bin_from+1, "interacting bins:",  len(int_strength))
	bins_sorted_by_strength = [b for b in sorted(int_strength, key=int_strength.get, reverse=True)]

	###############################
	# regions we are interested in (regions that contain TF binding sites)
	ranges = read_ranges(tfbs_file, tf_name)

	###############################
	# output: for each bin that interacts with promoter,
	# print all tfbs (actually ChIPSeq) region that fall within that bin
	hic_file = "{}/{}_{}_{}.tsv".format(hic_dir, gene_name, tf_name, assembly)
	outf = open(hic_file,"w")
	outf.write("\t".join(["% name", "chipseq_region", "hic_region", "hic_int_score"])+ "\n")
	rank = 0
	for b in bins_sorted_by_strength:
		bfrom, bto =bin_positions[b][1:3]
		if tad_start<=bfrom<=tad_end  or  tad_start<=bto<=tad_end:
			rank += 1
			tfbs_inside = find_tfbs_in(bfrom, bto, ranges)
			for t in tfbs_inside:
				outf.write( "\t".join( [tf_name, t, "{}_{}".format(bfrom, bto), "%5.0f"%(int_strength[b])]) + "\n")


	outf.close()

	return True


#########################################
########################################
if __name__ == '__main__':
	main()
