#!/usr/bin/python3

#
# This file is part of Progesterone pipeline.
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
from utils.utils import *
from utils.mysqldb import *



# this will store only bins within the TAD that contains our gene of interest
#########################################
#
def main():

	assembly  = "hg19"
	gene_name = "Hand2"
	species = "human"
	external_exp_id = "ENCFF633ORE"
	conf_file = "/home/ivana/.mysql_conf"

	contacts_file = "/storage/databases/encode/ENCSR551IPY/ENCFF331ABN.h5"

	for f in [contacts_file,conf_file]:
		if not os.path.exists(f):
			print(f,"not found")
			exit()
	db = connect_to_mysql(conf_file)
	cursor = db.cursor()
	switch_to_db(cursor,'progesterone')
	search_db(cursor,"set autocommit=1")

	# find xref_id for the experimental data file
	exp_file_xref_id = get_xref_id(db,cursor,external_exp_id)

	# find gene coordinates
	[chromosome, strand, gene_start, gene_end] = get_gene_coords (db,cursor, gene_name, assembly)

	# ind the TAD
	[tad_start, tad_end] = get_tad_region(db, cursor, exp_file_xref_id, chromosome, gene_start, gene_end)


	###############################
	# for promoter range take 1Kbp from the tx start region
	# it doesn't really matter, because the resolution here
	# is determined by HiC resolution, which is 40 Kbp
	promoter_range = []
	if strand=='+':
		promoter_range = [gene_start-1000, gene_end]
	elif strand=='-':
		promoter_range = [gene_start, gene_end+1000]
	else:
		print ("strand %s not recognized"%strand)
		exit()

	print("%s  chromosome: %s   promoter range:" % (gene_name, chromosome), promoter_range)


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
	tad_bins = []
	for b in range(chrom_bin_from,chrom_bin_to+1):
		bfrom, bto =bin_positions[b][1:3]
		if bfrom<=pstart<=bto or  bfrom<=pend<=bto: promoter_bins.append(b)
		if tad_start<=bfrom<=tad_end or tad_start<=bto<=tad_end: tad_bins.append(b)

	# generalize, at some point
	if len(promoter_bins)>1:
		print("bins containining promoter:", promoter_bins)
		for pb in promoter_bins:
			print (pb, infile['bin_positions'][pb])
		print ("promoter straddles two intervals; please generalize me accordingly")
		exit()
	pb = promoter_bins[0]
	print("promoter bin: ", pb)
	print("numer of tad bins", len(tad_bins))

	###############################
	# interactions between the promoter-containing bin and everybody else
	pb_interactions = infile['interactions'][pb][:]
	int_strength = {}
	for b in  tad_bins:
		if not math.isnan(pb_interactions[b]) and pb_interactions[b]>0.0:
			int_strength[b] = pb_interactions[b]

	###############################
	# store hic region that contains  our gene
	bfrom, bto = bin_positions[pb][1:3]
	fields  = {'species':species, 'chromosome':chromosome, 'assembly':assembly, 'rtype':'interacting',
						'rfrom':bfrom, 'rto':bto, 'xref_id':exp_file_xref_id}
	gene_hic_region_id = store_or_update(cursor, 'regions', fields, None)

	for b in  tad_bins:
		bfrom, bto =bin_positions[b][1:3]
		# store region (address)
		if b==pb:
			interacting_hic_region_id = gene_hic_region_id
		else: # we've done it already
			fields  = {'species':species, 'chromosome':chromosome, 'assembly':assembly, 'rtype':'interacting',
							'rfrom':bfrom, 'rto':bto, 'xref_id':exp_file_xref_id}
			interacting_hic_region_id = store_or_update(cursor, 'regions', fields, None)

		# store info about the binding region
		fields = {'gene_name':gene_name, 'gene_hic_region_id':gene_hic_region_id,
				  'interacting_hic_region_id':interacting_hic_region_id, 'interaction':int_strength[b] }
		interaction_id = store_or_update (cursor, 'hic_interactions', fields, None)



	###############################
	return True


#########################################
########################################
if __name__ == '__main__':
	main()
