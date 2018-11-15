#!/usr/bin/python3


# use Biomart to get predictiond for Hand2 only;
# careful -- need hg19 here http://grch37.ensembl.org/biomart
# range: 173880001, 175320000 (from 02_emve_tads.py for Hand2)

# to convert between assemblies: (CrossMap.py in python utils) for example
# CrossMap.py bed /storage/databases/liftover/hg38ToHg19.over.chain test.bed3
# bed3 consists of 3 cols only (chrom from to) but the first entry must be a string (chr1, chr2 etc)
# CrossMap website with chainfiles:  http://crossmap.sourceforge.net/#chain-file

