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
from statistics import mean, stdev
from math import floor
import matplotlib.pyplot as plt


def rescaled(x, offset, scale):
    return (x - offset) / scale

#########################################
def plot_1(ax, starts, ends, offset, scale):
    barwidth = 20
    y_offset = -barwidth
    prev_level_0_end = -1
    for i in range(len(starts)):
        [start, end] = [rescaled(starts[i], offset, scale), rescaled(ends[i], offset, scale)]
        if start>prev_level_0_end:
            prev_level_0_end = end
            y_offset = -barwidth
        y_offset += barwidth
        ax.hlines(y= y_offset, xmin=start, xmax=end, linewidth=barwidth*3, color='b')

#########################################
def plot_2(ax, starts, ends, offset, scale):
    barwidth = 20

    number_of_bins = 1000
    bin_counts = [0]*number_of_bins
    bin_positions = [0.0]*(number_of_bins)
    bin_positions[0] = 1.0/number_of_bins/2
    for b in range(1,number_of_bins): bin_positions[b] = bin_positions[b-1]+1.0/number_of_bins

    for i in range(len(starts)):
        [start, end] = [rescaled(starts[i], offset, scale), rescaled(ends[i], offset, scale)]
        first_bin = floor(start*number_of_bins)
        last_bin  = floor(end*number_of_bins)
        for b in range(first_bin,last_bin):
            bin_counts[b]+=1

    ret = ax.bar(bin_positions, bin_counts, width=1.0/number_of_bins)
    print(ret)

#########################################
def main():
    inpath = "raw_data/tads"
    chromosome = 1

    infile = "{}/chr{}.tsv".format(inpath, chromosome)
    for d in [inpath, infile]:
        if not os.path.exists(d):
            print(d, "not found")
            exit()

    starts = []
    ends = []
    s_stdvs = []
    e_stdvs = []
    inf = open(infile, "r")
    for line in inf:
        if line[0] == '%': continue
        [start, stdv_start, end, stdv_end] = line.split("\t")[1:5]
        starts.append(int(start))
        s_stdvs.append(int(stdv_start))
        ends.append(int(end))
        e_stdvs.append(int(stdv_end))
    inf.close()

    offset = min(starts)
    scale = max(ends) - offset
    fig, ax = plt.subplots(nrows=2,ncols=1)

    plot_1 (ax[0], starts, ends, offset, scale);
    plot_2 (ax[1], starts, ends, offset, scale);
    plt.show()

    return

#########################################
if __name__ == '__main__':
    main()
