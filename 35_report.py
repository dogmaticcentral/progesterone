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
	return chipseq_regions


#########################################
def get_tfbs_files(binding_site_dirs, assembly, gene_name):
	fnms = []
	for tfbs_dir in binding_site_dirs:
		for fnm in os.listdir(tfbs_dir):
			if not fnm[-4]!=".tsv": continue
			if not assembly in fnm: continue
			if not gene_name in fnm: continue
			fnms.append("{}/{}".format(tfbs_dir,fnm))
	return fnms

def get_int_strength(hic_dirs, assembly, gene_name, tf_name, region):
	intstr = ""
	interaction = {}
	for hic_dir in hic_dirs:
		for fnm in os.listdir(hic_dir):
			if not fnm[-4]!=".tsv": continue
			if not assembly in fnm: continue
			if not gene_name in fnm: continue
			if tf_name in fnm:
				intfor = tf_name
				cmd = "grep {} {}/{}".format(region, hic_dir, fnm)
			elif "self" in fnm:
				intfor = "self"
				cmd = "grep {} {}/{}".format("self", hic_dir, fnm)
			else:
				continue
			for line in subprocess.getoutput(cmd).split("\n"):
				print(cmd)
				print(line)
				interaction[intfor] = int(line.rstrip().split().pop())
				# if there is more than one line ... too bad
				break

	if tf_name in interaction or  "self" in interaction:
		intstr = "%d%%"%(int(interaction[tf_name]/interaction["self"]*100))
	return intstr

#########################################
def process(dirs,  assembly, gene_name, tf_name, region):
	# if interaction strength available, spit it out
	int_strength = get_int_strength(dirs["hic"],  assembly, gene_name, tf_name, region)
	print("\t{}  {}".format(region, int_strength))


#########################################
def report (dirs, assembly, chromosome, gene_name, tf_name):

	for fnm_path in get_tfbs_files(dirs["binding_site"], assembly, gene_name):
		chipseq_regions = read_tfbs_ranges(fnm_path, chromosome,  tf_name)
		if len(chipseq_regions)==0: continue
		print(fnm_path)
		for region in chipseq_regions:
			process(dirs, assembly, gene_name, tf_name, region)

#########################################
def main():

	assembly   = {'human':'hg19', 'mouse':'mm9'}
	gene_name  = "Hand2"
	chromosome = {'human':'4', 'mouse':'8'}
	dirs = {}
	dirs["binding_site"] = ['raw_data/tf_binding_sites_geo', 'raw_data/tf_binding_sites_ucsc']
	dirs["hic"] = ['raw_data/hic_interactions']
	for species in ['human','mouse']:
		for tf_name in ['ESR1', 'PGR']:
			report (dirs, assembly[species], chromosome[species], gene_name, tf_name)



#########################################
########################################
if __name__ == '__main__':
	main()
