#!/usr/bin/python3

import subprocess, os, urllib.request
from bs4 import BeautifulSoup
from Bio import motifs
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
import numpy as np


#########################################
def read_bed(infile, region_chrom, region_start, region_end):
	intervals = []
	if region_start: region_start = int(region_start)
	if region_end:   region_end   = int(region_end)
	inf = open(infile, "r")
	for line in inf:
		fields = line.rstrip().split("\t")
		if len(fields)<3: continue
		chrom = fields[0].replace('chr','')
		if chrom != region_chrom: continue
		try:
			[start,end] = [int(i) for i in fields[1:3]]
		except:
			continue
		# pass None as region start/end to skip this filter
		if (region_start and region_end) and (end<=region_start or region_end<=start): continue
		intervals.append([start,end])
	inf.close()
	return intervals


#############
def overlap(interval_list, qry_interval):
	ovlp = False
	for interval in interval_list:
		if qry_interval[1]<interval[0]: continue
		if interval[1]<qry_interval[0]: continue
		# some leeway could be left here one day ...
		ovlp = True
	return ovlp


#############
def read_binding_intervals(data_dir, agonist_file, vehicle_file, chrom, region_start, region_end):
	# agonist
	infile = "{}/{}".format(data_dir, agonist_file)
	agonist_binding_intervals = read_bed(infile, chrom, region_start, region_end)

	if not vehicle_file: return agonist_binding_intervals

	# if we have control file, subtract regions that popped up
	# with vehicle only:
	infile = "{}/{}".format(data_dir, vehicle_file)
	vehicle_binding_intervals = read_bed(infile, chrom, region_start, region_end)

	for interval in agonist_binding_intervals:
		if overlap(vehicle_binding_intervals, interval):
			agonist_binding_intervals.remove(interval)
			continue

	return agonist_binding_intervals


#########################################
def read_pfm(jaspar_motifs_file, tf_name):
	motif = None
	with open(jaspar_motifs_file) as handle:
		for m in motifs.parse(handle, "jaspar"):
			if m.name == tf_name:
				motif = m
				break
	return motif


#########################################
def read_or_download_sequence(chipseq_regions_dir, assembly, chromosome, tf, start, end):
	seqfile = "{}/{}_{}_{}_{}_{}.txt".format(chipseq_regions_dir, tf, assembly, chromosome, start, end)
	if (os.path.exists(seqfile) and os.path.getsize(seqfile) > 0):
		outf = open(seqfile, "r")
		seq = outf.read()
		outf.close()
	else:
		seq = ucsc_fragment_sequence(assembly, chromosome, start, end)
		outf = open(seqfile, "w")
		outf.write(seq.replace("\n", ""))
		outf.close()

	return seq



#########################################
def ucsc_fragment_sequence(assembly, chrom, start, end):
	if not 'chr' in chrom: chrom = 'chr'+chrom
	das_request = "http://genome.ucsc.edu/cgi-bin/das/{}/".format(assembly)
	das_request += "dna?segment={}:{},{}".format(chrom, start, end)
	response = urllib.request.urlopen(das_request)
	html = response.read()
	soup = BeautifulSoup(html, 'html.parser')
	if not soup: return None
	return soup.find('dna').string.strip().replace("\n", "")


#########################################
def ucsc_gene_coords(gene_name, ucsc_gene_regions_dir):
	cmd = "grep -i %s %s/*" % (gene_name, ucsc_gene_regions_dir)
	ret = subprocess.check_output(cmd, shell=True).decode('utf-8').rstrip()
	if not ret or len(ret) == 0:
		print("no entry for %s found in %s " % (gene_name, ucsc_gene_regions_dir))
		exit()

	lines = []
	for line in ret.split("\n"):
		fields = line.split("\t")
		[infile, refseq_id] = fields[0].split(":")
		lines.append(fields[1:])
	if len(lines) == 0:
		print("no entry for %s found in %s " % (gene_name, ucsc_gene_regions_dir))
		return None, None
	if len(lines) == 2:
		print("more than one entry found for %s found in %s " % (gene_name, ucsc_gene_regions_dir))
		return None, None
	# we assume certain format in the file name, containing the chromosome number: e.g. chr18.csv
	chromosome = infile.split("/")[-1].replace(".csv", "")
	[name, strand, txStart, txEnd] = lines[0]
	return chromosome, strand, [int(txStart), int(txEnd)]

#########################################
def parsed_alignment(maf_file):
	infile = open(maf_file, "r")
	labels = []
	seqs = {}
	for line in infile:
		if line[0] != 's': continue
		fields = line.rstrip().split()
		if len(fields)!=7: continue
		[src, start,size,strand, src_size, sequence] = fields[1:7]
		[start,size,src_size] = [int(i) for i in [start,size,src_size]]
		fields = src.split(".")
		assembly = fields[0]
		chrom = ".".join(fields[1:]).replace("_","")
		if strand=='+':
			rfrom = start
			rto = start + size -1
		else:
			rto   = src_size - start + 1
			rfrom = rto - size + 1
		label = "_".join([assembly, chrom, str(rfrom), str(rto), strand])
		if not label in labels:
			labels.append(label)
			seqs[label] = ""
		seqs[label] += sequence
	infile.close()

	return labels, seqs

#########################################
def get_alignment(species, assembly, chrom, region_from, region_to, scratch):
	# from http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_6
	maf_region_extraction_tool = "/usr/bin/mafsInRegion"
	# mafs come from here http://hgdownload.cse.ucsc.edu/downloads.html
	# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/multiz100way/
	# or whichever species or assembly appropriate
	if not 'chr' in chrom: chrom = "chr"+chrom
	maf_file = "/storage/databases/ucsc/mafs/{}/{}/{}.maf".format(species, assembly, chrom)

	for dep in [maf_region_extraction_tool, maf_file, scratch]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()
	if not os.path.isdir(scratch):
		print(dir, "is not directory")
		exit()

	bed_in  = "{}/{}.bed".format(scratch, os.getpid())
	maf_out = "{}/{}.maf".format(scratch, os.getpid())
	with open(bed_in,"w") as outf:
		# it looks like  mafsInRegion counts frmo 0 and takes the region_to as a non-inclusive upper limit
		outf.write("{} {} {}\n".format(chrom, region_from-1, region_to))

	cmd = "{} {} {} {}".format(maf_region_extraction_tool, bed_in, maf_out, maf_file)
	subprocess.call(cmd, shell=True)

	labels, sequences = parsed_alignment(maf_out)

	os.remove(maf_out)
	os.remove(bed_in)

	return labels, sequences


#########################################
def remove_all_gaps(almt):
	almt_length = min([len(seq) for seq in almt.values()])
	gapped_pos = [i for i in range(almt_length) if len([seq for seq in almt.values() if seq[i]!='-'])==0]
	for nm,seq in almt.items():
		almt[nm] = ''.join([seq[i] for i in range(almt_length) if not i in gapped_pos ])


#########################################

#########################################
if __name__ == "__main__":
	# print(ucsc_fragment_sequence('mm10',8, 57805369, 57805386))
	#print(ucsc_gene_coords('Hand2', "/storage/databases/ucsc/gene_ranges/mouse/mm9"))
	##############################
	# careful parsing maf:
	# from https://genome.ucsc.edu/FAQ/FAQformat.html#format5
	# The "s" lines together with the "a" lines define a multiple alignment; the fields
	# are src | start | size | strand | srcSize | sequence
	# The start of the aligning region in the source sequence  is a zero-based number.
	# If the strand field is "-" then this is the start relative to the reverse-complemented source sequence.
	# That is - it should be subtracted from the size column
	# https://raw.githubusercontent.com/bxlab/bx-python/master/scripts/maf_to_fasta.py
	# does not take this into account

	# /usr/bin/maf_parse -s 58814511 -e 58814525 /storage/databases/ucsc/mafs/mouse/mm9/chr8.maf
	# a score=938862.000000
	# s mm9.chr8                        58814510   15 + 131738871 G-AT-GC-ATTTTGTCTT
	# s cavPor2.scaffold_280174            91854   13 +     98733 A-AT-GA-ATTACATC--
	# q cavPor2.scaffold_280174                                   9-99-99-99999999--
	# s rn4.chr16                       52832882   15 -  90238779 G-AT-AC-CTTATGTCTT
	# q rn4.chr16                                                 9-99-99-9999999999
	# s hg18.chr4                       15657414   15 - 191273063 T-AT-GA-ATTCTGTCTT
	# s panTro2.chr4                    15951551   15 - 194897272 T-AT-GA-ACTCTGTCTT

	# incidentially we discover that mafsInRegion (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/mafsInRegion)
	# is significantly faster

	get_alignment('mouse', 'mm10', 'chr8', 57805369, 57805386, '/home/ivana/scratch', '/home/ivana/scratch/test.afa')




