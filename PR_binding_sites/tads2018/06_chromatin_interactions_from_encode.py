#!/usr/bin/python3

# source https://www.encodeproject.org/experiments/ENCSR551IPY/
# hdf files can be inspected with h5dump --contents <filename>
# even easier: h5ls <filename>  (or h5ls -vlr <filename> but this might be too verbose)
import h5py

def read_tfb_ranges(infile):
	ranges = {}
	chrom = None
	for line in open(infile, "r"):
		[encodebin, chrom_readin, chromStart, chromEnd, name, score, expCount, expNums, expSCores] = line.rstrip().split("\t")
		if chrom==None:
			chrom=chrom_readin
		elif (chrom != chrom_readin):
			print ("Unexpected: chrom not the same for all entry lines in %s", infile)
			exit()
		if not name in ranges: ranges[name]=[]
		ranges[name].append([int(chromStart), int(chromEnd)])
	return ranges

#########################################
def main():

	# regions we are interested in (expected usage: from the previous script, 04_tf_binding_sites_from_UCSC.py)
	infile = "raw_data/Hand2_tfbs.tsv"

	ranges = read_tfb_ranges(infile)

	#print(ranges['ESR1'])

	contact_file= "/storage/databases/encode/ENCSR551IPY/ENCFF331ABN.h5"
	infile= h5py.File(contact_file,'r')


	return True


#########################################
########################################
if __name__ == '__main__':
	main()
