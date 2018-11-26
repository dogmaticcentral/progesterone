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

# hdf files can be inspected with h5dump --contents <filename>
# even easier: h5ls <filename>  (or h5ls -vlr <filename> but this might be too verbose)
import h5py

#########################################
def main():


	###################
	# small  file
	contact_file= "/storage/databases/encode/ENCSR551IPY/ENCFF354SCA.h5"
	infile= h5py.File(contact_file,'r')

	print("\ndatasets in %s"%contact_file)
	for k in infile.keys():
		print (k)

	# the info in 'chromosomes' is the list of chromosome names.
	# 'chr1'  with index of 0, through 'chrM' with index of 24
	print("\nchromosomes")
	print (infile['chrs'])
	print (infile['chrs'].dtype)
	print (infile['chrs'][...])
	print (infile['chrs'][1:3])

	# 'bin_positions' is a list of bin descriptors, giving [chrom, from, to]
	# for each interval (bin) for which the interaction will be presented in 'interactions'
	print("\nbin_positions")
	print (infile['bin_positions'])
	print (infile['bin_positions'][320:][:])

	# 'chromosome ranges' in principle provides the same info as in 'bin_positions'
	# but organized differently (this is not SQL, queries have to be anticipated)
	# it gives, for each chromosome, the [bin_from, bin_to] range in the bin_positions list
	print("\nchromosome ranges")
	print (infile['chr_bin_range'])
	print (infile['chr_bin_range'][:5][:])

	# find the interactions within chromosome Y:
	print("\nchromosome Y interactions")
	chromY_index = 23
	print ("chrom Y index and name: ", chromY_index, infile['chrs'][chromY_index])
	print ("chrom Y bins: ", infile['chr_bin_range'][chromY_index])
	[bin_from, bin_to] = infile['chr_bin_range'][chromY_index]
	# not that bin_to is inclusive, not the upper limit
	for bin_pos in infile['bin_positions'][bin_from:bin_to+2]:
		print(bin_pos)
	print('interactions')
	for idx1 in range(bin_from, bin_to+1):
		print (infile['interactions'][idx1][bin_from:bin_to+1])

	# note we can also pull out the block
	interaction_block = infile['interactions'][bin_from:bin_to+1][bin_from:bin_to+1]

	infile.close()

	###################
	# big  file also contains
	contact_file= "/storage/databases/encode/ENCSR551IPY/ENCFF331ABN.h5"
	infile= h5py.File(contact_file,'r')

	# this one also contains balance_factors
	print("\n\ndatasets in %s"%contact_file)
	for k in infile.keys():
		print (k)
	# how long woudl it take to pull out the interaction block now?
	# about 1 sec for chrom Y (chrom_index = 23)
	# about 5 sec for chrom 4 (chrom_index =  3)
	# about 6.5 sec for chrom 1 (chrom_index =  0)
	chrom_index = 3
	[bin_from, bin_to] = infile['chr_bin_range'][chrom_index]

	interaction_block = infile['interactions'][bin_from:bin_to+1][bin_from:bin_to+1]

	infile.close()



	return True


#########################################
########################################
if __name__ == '__main__':
	main()
