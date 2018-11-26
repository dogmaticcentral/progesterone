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

# ita appears that the interactoin matrices in ENCSR551IPY are already balanced,
# in the sense advertised in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4347522/
# "This procedure attempts to balance the matrix by equalizing the sum of every row/column in the matrix."

import h5py
import math, sys
#########################################
def main():

	if len(sys.argv) < 2:
		print  ("usage: %s <'rows'|'cols'|'both'>" % sys.argv[0])
		exit()
	check_both = False
	check_rows = True
	if sys.argv[1]=='rows':
		check_rows = True
	elif sys.argv[1]=='cols':
		check_rows = False
	elif sys.argv[1]=='both':
		check_both = True
	else:
		print("unknown argument: %s" % sys.argv[1])
		exit()

	###################
	# big  file also contains
	#contact_file= "/storage/databases/encode/ENCSR551IPY/ENCFF354SCA.h5" # smallest - not normalized
	#contact_file= "/storage/databases/encode/ENCSR551IPY/ENCFF331ABN.h5" #  the biggest
	contact_file= "/storage/databases/encode/ENCSR551IPY/ENCFF383OVH.h5" # second smallest
	#contact_file= "/storage/databases/encode/ENCSR551IPY/ENCFF430DYJ.h5" # second biggest

	print ("checking", contact_file)

	infile= h5py.File(contact_file,'r')
	#infile['bin_positions'].shape is a tuple (nr_bins, pieces_of_data_per_bin)
	# pieces_of_data_per_bin is always 3: the bin address [chrom, from, to]
	number_of_bins = infile['bin_positions'].shape[0]

	# row sums before normalization
	# caveat - some entries are nana
	ints =  infile['interactions'] # this is hdf4 dataset 0 it can be slced like  ints[10,:] or  ints[:,10]
	
	if check_rows or check_both:
		row_sum_0 = None
		for i in range(number_of_bins):
			if i%100==0: print("row",i,"out of",number_of_bins )
			row_sum = sum ( [x  for x in ints[i,:] if not math.isnan(x) and x>0.0 ] )
			if not row_sum>0 : continue
			if not row_sum_0: row_sum_0= row_sum
			if  abs(row_sum_0-row_sum)>1:
				print (i, row_sum)
				#print ([x  for x in ints[i] if not math.isnan(x) and x>0.0])
				#exit()
		print("row sum:",row_sum_0)

	if (not check_rows) or check_both:
		col_sum_0 = None
		for j in range(number_of_bins):
			if j%100==0: print("col",j,"out of",number_of_bins )
			col_sum = sum ( [x  for x in ints[:,j] if not math.isnan(x) and x>0.0 ] )
			if not col_sum>0 : continue
			if not col_sum_0: col_sum_0= col_sum
			if  abs(col_sum_0-col_sum)>1:
				print (j, col_sum)
	if check_both:
		print("row sum:",row_sum_0)
		print("col sum:",col_sum_0)


	infile.close()



	return True


#########################################
########################################
if __name__ == '__main__':
	main()
