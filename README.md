# Progesterone

This is a set of scritps that try to answer the question: is it possible that a gene A 
is under control of transcription  factor B, given currently available experimental 
evidence. Response to progesterone and estrogen was the original topic.

_Progesterone_ is not a library  - the scripts are loosely connected by a couple of common methods 
in the utils module. They average some hundred lines in length, and can
probably be figured out without much explanation. 
They have been tested to run under python 3.6.6, and are supposed 
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
  
* ChIPSeq regions from ENCODE experiment, [collected in UCSC genome db](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeRegTfbsClustered/)

* [CrossMap](http://crossmap.sourceforge.net/) for transforming coordinates to a single reference assembly 
  (together with transformation chain files, [here](http://crossmap.sourceforge.net/#chain-file))

* [Gnuplot](http://www.gnuplot.info/)

* Matplotlib (_pip3 install matplotlib_; 
if your scripts complain about missing tkinter, you might also have to do
_sudo apt install python3-tk_)

* [tools](https://www.h5py.org/) for handling data in [HDF5 format](https://portal.hdfgroup.org/display/support)

* TF binding motif databases [JASPAR](http://jaspar.genereg.net/) and 
   [Hocomoco](http://hocomoco11.autosome.ru/)

* [Motif](http://biopython.org/DIST/docs/api/Bio.motifs-module.html) 
  from [Biopython](https://biopython.org/)

* more ChIPSeq data from [GEO](https://www.ncbi.nlm.nih.gov/geo/) database ([GSE34927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34927),
 [GSE36455](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE36455), for example. Again YMMV.)

## Gene coordinates
You may start by downloading chromosome lengths and gene coordinates from UCSC using 
[01_chromosome_lengths_from_UCSC.py](01_chromosome_lengths_from_UCSC.py) 
and [02_gene_ranges_from_UCSC.py](02_gene_ranges_from_UCSC.py), 
or you can download them later, when they become needed.

## What's with this TAD business


With 2 meters of DNA squished  in the nucleus with 6 micrometres in diameter, we expect 
that the packing is quite complex. What we know in addition, is that 
the packing is not random. Rather, the genetic material is organized into domains, 
[TADs](https://en.wikipedia.org/wiki/Topologically_associating_domain),
that stay localized in space, even as they  keep rearranging internally, and switching their positions within the nucleus, reviewed in 
[Merkenschlager & Nora](https://www.researchgate.net/profile/Elphege_Nora/publication/301482856_CTCF_and_Cohesin_in_Genome_Folding_and_Transcriptional_Gene_Regulation/links/5726c5b508ae262228b21511/CTCF-and-Cohesin-in-Genome-Folding-and-Transcriptional-Gene-Regulation.pdf).
The implication would be that any transcription factor (TF) binding site therein is potentially regulating
any (or all) genes therein. This means that we can look for the regulating TF binding sites as far as ~MBp away from 
the promoter, an idea that would be an anathema half a decade back. In addition it opens a happy possiblity that genes 
are not under a single TF site control, but under control of a number of them, in a stochastic way, with
the probablity of TF having an effect being proportional to the number of TF binding sites in a TAD.

Where on the chromosome, then, is the TAD that my gene belongs to?
The idea repeatedly appears in the literature that TAD boundaries are conserved across cell types and species
(for example, [Dixon _et al_ review, Moll Cell 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5371509/)). 
That would be useful, beacuse with unversal definetion of TADs  we could write very many scripts with very many 
purposes, having looked up the TAD definition only once.

 [10_tads_overview.py](10_tads_overview.py) explores that possibility  - with mixed results.

You will need to download TAD files from the [Yue lab page](http://promoter.bx.psu.edu/hi-c/publications.html) 
 (we suggest  sticking with hg19 throughout the pipeline  - use link named 'TADs in  hg19'). 
 Adjust the dirpath in  [10_tad_overview.py](10_tad_overview.py) accordingly. This script groups TAD intervals 
 that appear repeatedly in different experiments, in the hope of finding those
 that are supposedly conserved across different cell types. It seems that one could make some progress
 looking for regions that are rarely assigned to a domain, and thus delineate the topmost division in a 
 TAD hierarchy. The number of these divisions, however, appears to be an order of magnitude smaller
 than the number of TADs typically reported in each experiment.
 
 [11_tads_pic.py](11_tads_pic.py) contains basic code to illustrate  the reported TAD domains. 
 It uses [Matplotlib](https://matplotlib.org/); note that in Matplotlib you can zoom 
 into any region of the graph you are interested in. 
 Here are the TADS for human chromsome 1, from 35 different experiments from the Yue lab collection:
 
![seepic](images/tad_exp_chr1.png)
 
 The whole chromosome length has been rescaled to the range [0,1].
If you squint a little you can see that the basic TAD  structure corresponds to the regions of densest
gene occupation region (middle and bottom panels) -  careful with the interpretation of this graph: see 
[11_tads_pic.py](11_tads_pic.py):plot_2()). However the individual TAD assignments 
(top panel; each horizontal level corresponds to one experiment) 
vary widely between different experiments. Therefore we choose to stick with a single experiment
([Homo sapiens endometrial microvascular endothelial cells](https://www.encodeproject.org/experiments/ENCSR551IPY/))
because the  cell type it uses  matches most closely the type of cells we are interested in.

## Which TAD does my gene belong to

## Where are ChIPSeq regions on the chromosome
_Progesterone_ pipeline takes hg19 and mm9  as its reference assemblies. 
In GEO repositories there is usually a file called *_series_matrix.txt 
(which is actually a tsv file), where this info can be found. 
You may grep for hg or mm to see which one is referred to. If it is not hg19 for human
or mm9 for mouse, the coordinates need to be translated.

The scripts expect "bed" format, which here means that the columns are tab separated, the
first column is chromosome number (possibly prefixed by 'chr'), and the following two are
region start and region end. If the file is not
too far from that format, you can perhaps help yourself out with 
[linux _cut_ command](https://www.thegeekstuff.com/2013/06/cut-command-examples/), 
see also 
[here](https://unix.stackexchange.com/questions/35369/how-to-define-tab-delimiter-with-cut-in-bash) 
for tab delimiter handling with _cut_. For example

`cut -d$'\t' -f3-5  GSM857545_1_PR_oil_s_4_aligned.tsv > GSM857545_1_PR_oil_s_4_aligned.bed`


 (In Windows you might try using a spreadsheet program to reformat the file. Just make sure
 you do not have too many  in your bed file.)

## Which regions come  in contact within the TAD (and how often)

## Deja vu all over again: mouse data

## Putting it all together

## So, functional TF sites or not?
Print the table you have just made, and take it to a competent experimentalist. 
(S)he might be able to help.
