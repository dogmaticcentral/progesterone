# Progesterone

This is a set of scritps that try to answer the question: is it possible that a gene A 
is under control of transcription  factor B, given currently available experimental 
evidence. Response to progesterone and estrogen was the original topic.

_Progesterone_ is not a library  - the scripts are loosely connected by a couple of common methods 
in the utils module. They have been tested to run under python 3.6.6, and are supposed 
to be run in the  order in which they are enumerated. They answer a series of 
sub-questions, which may add up to an answer. Depends on the answer you are hoping for.

## Dependencies

_Progesterone_ pipeline depends on several data sources and python packages. However, you can 
download/install them only when they become necessary. 
* [collection of TADs from Yue lab](http://promoter.bx.psu.edu/hi-c/publications.html)

* [UCSC genome database, accessed directly through MySQL](https://genome.ucsc.edu/goldenpath/help/mysql.html)

* MySQLdb, installed with _sudo apt install python3-mysqldb_

* data from HiC experiment on [... something](https://www.encodeproject.org/experiments/ENCSR551IPY/) - 
  adapt this to the cell/tissue type you are interested in
  
* ChIPSeq regions from ENCODE experriment, [collected in UCSC genome db](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/)

* [Gnuplot](http://www.gnuplot.info/)

* [tools](https://www.h5py.org/) for handling data in [HDF5 format](https://portal.hdfgroup.org/display/support)

* TF binding motif databases [JASPAR](http://jaspar.genereg.net/) and 
   [Hocomoco](http://hocomoco11.autosome.ru/)

* [Motif](http://biopython.org/DIST/docs/api/Bio.motifs-module.html) 
  from [Biopython](https://biopython.org/)

* more ChIPSeq data from [GEO](https://www.ncbi.nlm.nih.gov/geo/) database ([GSE34927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34927),
 [GSE36455](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36455), for example. Again YMMV.)


## What's this TAD business

## Which TAD does my gene belong to

## Where are ChIPSeq regions on the chromosome

## Which regions come  in contact within the TAD (and how often)

## Deja vu all over again: mouse data

## Putting it all together

## So, functional TF sites or not?
Print the table you have just made, and take it to a competent experimentalist. 
(S)he might be able to help.