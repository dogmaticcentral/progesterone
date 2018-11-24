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

# the input files with tad domains are from
# http://promoter.bx.psu.edu/hi-c/publications.html
# conclusion: very much against the claims of Dixon et al,
# the TADs are significantly different from
# one cell type to the next

import os
from math import floor
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection


#########################################
def rescaled(x, offset, scale):
	return (x - offset) / scale


#########################################
# plots input intervals from one file as strips on the same y-height in diagram
def plot_0(ax, dirpath, datafiles, chromosome, offset, scale):

	patches = []
	height = .5
	y_offset = 0
	number_of_files = len(datafiles)
	for file in datafiles:
		inf = open("{}/{}".format(dirpath,file),"r")
		y_offset += 1
		for line in inf:
			[chr, s, e] = line.rstrip().split()
			if chr!= "chr{}".format(chromosome): continue
			[start,end] = [rescaled(int(s)/1000, offset, scale), rescaled(int(e)/1000, offset, scale)]
			rect = mpatches.Rectangle(xy=(start, y_offset-height/2), width=(end-start), height=height)
			patches.append(rect)
		inf.close()

	ax.set_ylim(0, number_of_files+1)
	ax.set_xlim(0,1)
	ax.add_collection(PatchCollection(patches))


#########################################
# illustrates clusters of intervals - the thickness of each strip
# corresponds to the number of experiments
# the y-position has no intrinsic meaning - the strip should not overlap
def plot_1(ax, counts, starts, ends, offset, scale):

	y_offset = -1
	prev_level_0_end = -1
	patches = []
	max_counts = max(counts)
	for i in range(len(starts)):
		[start, end] = [rescaled(starts[i], offset, scale), rescaled(ends[i], offset, scale)]
		if start>prev_level_0_end:
			prev_level_0_end = end
			y_offset = -1
		y_offset += 1
		height = counts[i]/max_counts
		rect = mpatches.Rectangle(xy=(start, y_offset-height/2), width=(end-start), height=height)
		patches.append(rect)

	ax.set_ylim(0,10)
	ax.set_xlim(0,1)
	ax.add_collection(PatchCollection(patches))

#########################################
# number of interval clusters to which  each bin corresponds
def plot_2(ax, starts, ends, offset, scale):
	number_of_bins = 1000
	bin_counts = [0]*number_of_bins
	bin_positions = [0.0]*(number_of_bins)
	bin_positions[0] = 1.0/number_of_bins/2
	for b in range(1,number_of_bins): bin_positions[b] = bin_positions[b-1]+1.0/number_of_bins

	for i in range(len(starts)):
		[start, end] = [rescaled(starts[i], offset, scale), rescaled(ends[i], offset, scale)]
		first_bin = int(floor(start*number_of_bins))
		last_bin  = int(floor(end*number_of_bins))
		for b in range(first_bin,last_bin):
			bin_counts[b]+=1

	ax.set_xlim(0,1)
	ax.bar(bin_positions, bin_counts, width=1.0/number_of_bins)

#########################################
def read_processed_data(infile):
	counts  = [] # number of intervals grouping into this cluster
	starts = []
	ends   = []
	s_stdvs = []
	e_stdvs = []
	inf = open(infile, "r")
	for line in inf:
		if line[0] == '%': continue
		[count, start, stdv_start, end, stdv_end] = line.split("\t")[:5]
		counts.append(int(count))
		starts.append(int(start))
		s_stdvs.append(int(stdv_start))
		ends.append(int(end))
		e_stdvs.append(int(stdv_end))
	inf.close()

	return counts, starts, ends, s_stdvs, e_stdvs


#########################################
def  get_chrom_length(chr_lengths_file, chromosome):
	length = None
	chrstring  = "chr{}".format(chromosome)
	with open(chr_lengths_file,"r") as inf:
		for line in inf:
			[chrom,size] = line.split("\t")
			if chrstring!=chrom: continue
			try:
				length = int(size.strip())
			except:
				print ("error parsing chrom length line:", line)
				exit()
			finally:
				break
	if length==None:
		print ("Length for {} not found in {}".format(chromosome, chr_lengths_file))
		exit()
	return length

#########################################
def main():

	assembly   = "hg19"
	chromosome = '4'

	dirpath = "/storage/databases/3dgenomebrowser/hg19_TADs"
	inpath = "raw_data/tads"
	chr_lengths_file = "/storage/databases/ucsc/chromosome_lengths/%s.tsv" % assembly
	infile = "{}/chr{}.tsv".format(inpath, chromosome)
	for d in [inpath, infile, chr_lengths_file]:
		if not os.path.exists(d):
			print(d, "not found")
			exit()
	[counts, starts, ends, s_stdvs, e_stdvs ]= read_processed_data(infile)

	chrom_length = get_chrom_length(chr_lengths_file, chromosome)

	datafiles = []
	for path, dirs, files in os.walk(dirpath):
		datafiles += [file for file in files if file[-3:] == "txt"]

	#offset = min(starts)
	#scale = max(ends) - offset

	fig, ax = plt.subplots(nrows=3,ncols=1)
	plot_0(ax[0], dirpath, datafiles, chromosome,  0, chrom_length/1000)
	plot_1(ax[1], counts, starts, ends, 0, chrom_length/1000);
	plot_2(ax[2], starts, ends, 0,  chrom_length/1000);
	plt.show()

	return

#########################################
if __name__ == '__main__':
	main()
