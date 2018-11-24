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


# the single input file from endometrial microvascular endothelial cells from
# Job Dekker lab, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105710
# the (bed) file with TADs can be found here
# https://www.encodeproject.org/experiments/ENCSR551IPY/
# (under processed data)

from utils.utils import *

#########################################
def main():

	gene_name = "Hand2"

	tadfile = "/storage/databases/tads/encode/ENCFF633ORE.bed"
	ucsc_gene_regions_dir = "/storage/databases/ucsc/gene_ranges/human/hg19"

	for prerequisite in [tadfile, ucsc_gene_regions_dir]:
		if os.path.exists(prerequisite): continue
		print(prerequisite, "not found")
		exit()

	chromosome, strand, gene_range = ucsc_gene_coords(gene_name, ucsc_gene_regions_dir)
	if not chromosome or not gene_range:
		exit()

	print ( "{} {} {}:{}-{}".format(gene_name, strand, chromosome, gene_range[0], gene_range[1]) )

	[start, end] = get_tad (tadfile, chromosome, gene_range)
	print ("TAD containing %s region: %s:%d-%d   length %d"%(gene_name, chromosome, start, end, end-start))


#########################################
if __name__ == '__main__':
	main()




