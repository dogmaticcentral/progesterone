#!/usr/bin/python3
# the input files with tad domains are from
# http://promoter.bx.psu.edu/hi-c/publications.html
# conclusion: very much against the claims of Dixon et al,
# the TADs are significantly different from
# one cell type to the next

import os
from statistics import mean, stdev
from math import floor, ceil

#########################################
def find_place(interval_clusters, new_start, new_end):

    interval_placed = False
    new_length = new_end-new_start
    for interval_cluster in interval_clusters:
        for [start,end] in interval_cluster:
            length = end-start
            ds = abs(start-new_start)
            de = abs(end-new_end)
            if ds/length<0.2 and de/length<0.2 and ds/new_length<0.2 and de/new_length<0.2:
                interval_cluster.append([new_start, new_end])
                interval_placed = True
                break
        if interval_placed: break

    if not interval_placed:
        interval_clusters.append([[new_start, new_end]])
#############
def interval_stat(interval_cluster):

    if len(interval_cluster)==1:
        interval = interval_cluster[0]
        return [interval[0], 0, interval[1], 0, interval[1]-interval[0], 0]

    starts  = [interval[0] for interval in interval_cluster]
    ends    = [interval[1] for interval in interval_cluster]
    lengths = [interval[1]-interval[0] for interval in interval_cluster]

    return [int(i) for i in [mean(starts), stdev(starts), mean(ends), stdev(ends), mean(lengths), stdev(lengths)]]


def clean_up(tads):

    return

#########################################

def main():
    dirpath = "/storage/databases/3dgenomebrowser/hg19_TADS"

    datafiles = []
    for path, dirs, files in os.walk(dirpath):
        datafiles += [ file for file in files if file[-3:] == "txt"]

    tads = {}

    for file in datafiles:
        inf = open("{}/{}".format(dirpath,file),"r")
        for line in inf:
            [chr, start, end] = line.rstrip().split()
            #if chr !=  'chr22': continue
            if not chr in tads: tads[chr]=[]
            find_place(tads[chr], int(int(start)/1000), int(int(end)/1000))


    for chr, interval_clusters in tads.items():
        #print ("******************")
        #print(chr)
        stats = []
        for interval_cluster in interval_clusters:
            if len(interval_cluster)<5: continue
            stats.append([len(interval_cluster)] + interval_stat(interval_cluster))
        outf = open ("tads/{}.tsv".format(chr),"w")
        for stat in sorted(stats, key=lambda s: s[1]):
            outf.write("\t".join([str(s) for s in stat])+"\n")
        outf.close()

#########################################
if __name__ == '__main__':
    main()




