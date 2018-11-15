#!/usr/bin/python3

import subprocess, os, urllib.request
from bs4 import BeautifulSoup


#########################################
def ucsc_fragment_sequence(assembly,chrom,start,end):
	das_request  = "http://genome.ucsc.edu/cgi-bin/das/{}/".format(assembly)
	das_request += "dna?segment=chr{}:{},{}".format(chrom,start,end)
	response = urllib.request.urlopen(das_request)
	html = response.read()
	soup = BeautifulSoup(html,'html.parser')
	if not soup: return None
	return soup.find('dna').string.strip()

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
def get_tad(tadfile, chromosome, gene_range):
	tads = {}
	inf = open(tadfile, "r")
	for line in inf:
		[chr, start, end] = line.rstrip().split()[:3]
		if not chr in tads: tads[chr] = []
		tads[chr].append([int(start), int(end)])

	gene_tads = []
	for start, end in tads[chromosome]:
		if start <= gene_range[0] <= end or start <= gene_range[1] <= end:
			gene_tads.append([start, end])

	if len(gene_tads) == 0:
		print("TAD not found for chromosome", chromosome, "gene range", gene_range)
		exit()

	if len(gene_tads) > 1:
		print("gene straddles two TADS; generalize the code if needed")
		exit()

	return gene_tads[0]


#########################################
def get_alignment(species, assembly, chrom, region_from, region_to, scratch, outfile):
	# sudo apt install phast.v1_5.x86_64.deb (contains util maf_parse)
	maf_region_extraction_tool = "/usr/bin/maf_parse"
	# mafs come from here http://hgdownload.cse.ucsc.edu/downloads.html
	# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/multiz100way/
	# or whichever species or assembly appropriate
	maf_file = "/storage/databases/ucsc/mafs/{}/{}/{}/maf".format(species, assembly, chrom)
	# pip3 install bx-python
	# then put https://raw.githubusercontent.com/bxlab/bx-python/master/scripts/maf_to_fasta.py
	# into /usr/local/bin
	# make sure to change the shebang to python3 and make executable
	maf2afa_tool = "/usr/local/bin/maf_to_fasta.py"

	for dep in [maf_region_extraction_tool, maf_file, maf2afa_tool, scratch]:
		if not os.path.exists(dep):
			print(dep, "not found")
			exit()
	if not os.path.isdir(scratch):
			print(dir, "is not directory")
			exit()

	tmp_out = "{}/{}.maf".format(scratch,os.getpid())
	cmd = "{} -s {} -e {}  {} > {}".format(maf_region_extraction_tool,
											region_from, region_to, maf_file, tmp_out)
	subprocess.call(cmd)

	cmd = "{} < {} > {}".format(maf2afa_tool, tmp_out, outfile)
	subprocess.call(cmd)


	os.remove(tmp_out)

#########################################
if __name__=="__main__":
	print(ucsc_fragment_sequence('mm9',8,59791026,59791040))