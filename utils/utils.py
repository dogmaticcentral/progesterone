#!/usr/bin/python3

import subprocess, os, urllib.request
from bs4 import BeautifulSoup
from Bio import motifs
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import Seq
import numpy as np


#########################################
def read_pfm(jaspar_motifs_file, tf_name):
    motif = None
    with open(jaspar_motifs_file) as handle:
        for m in motifs.parse(handle, "jaspar"):
            if m.name == tf_name:
                motif = m
                break
    return motif


#########################################
def read_or_download_sequence(chipseq_regions_dir, assembly, chromosome, tf, start, end):
    seqfile = "{}/{}_{}_{}_{}_{}.txt".format(chipseq_regions_dir, tf, assembly, chromosome, start, end)
    if (os.path.exists(seqfile) and os.path.getsize(seqfile) > 0):
        outf = open(seqfile, "r")
        seq = outf.read()
        outf.close()
    else:
        seq = ucsc_fragment_sequence(assembly, chromosome, start, end)
        outf = open(seqfile, "w")
        outf.write(seq.replace("\n", ""))
        outf.close()

    return seq


#########################################
def almt_simplified(almtfile, pssm, revstrand):
    outf = open(almtfile, "r")
    assmbs = []
    seq = {}
    ranges = {}
    for line in outf:
        line = line.rstrip()
        if line[0] == '>':
            [asm, range] = line[1:].split(':')
            if not asm in seq:
                assmbs.append(asm)
                seq[asm] = ""
                ranges[asm] = []
            ranges[asm].append(range)
        else:
            seq[asm] += line
    outf.close()

    almt = ""
    for asm in assmbs:
        fields = asm.split(".")
        [species, chrom] = fields[:2]
        seq_straight = seq[asm].replace("-", "")[:20].upper()
        biopythonseq = Seq(seq_straight, unambiguous_dna)
        if revstrand:
            biopythonseq = biopythonseq.reverse_complement()
        try:
            score = pssm.calculate(biopythonseq)
            maxscore = np.amax(score)
        except:
            maxscore = -100
        almt += "%-10s %-20s %5.1f   %-15s %s\n" % (species, seq_straight, maxscore, chrom, ranges[asm])
        if (species == 'rn5'):
            almt += "-------------------------------------------\n"
    return almt


#####
def get_alignment_file(alignments_dir, species, assembly, chromosome, tf, start, end):
    almtfile = "{}/{}_{}_{}_{}_{}.txt".format(alignments_dir, tf, assembly, chromosome, start, end)
    if not (os.path.exists(almtfile) and os.path.getsize(almtfile) > 0):
        get_alignment(species, assembly, chromosome, start, end, alignments_dir, almtfile)
    return almtfile


#########################################
def ucsc_fragment_sequence(assembly, chrom, start, end):
    das_request = "http://genome.ucsc.edu/cgi-bin/das/{}/".format(assembly)
    das_request += "dna?segment=chr{}:{},{}".format(chrom, start, end)
    response = urllib.request.urlopen(das_request)
    html = response.read()
    soup = BeautifulSoup(html, 'html.parser')
    if not soup: return None
    return soup.find('dna').string.strip().replace("\n", "")


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

    if not gene_range:
        return tads[chromosome]

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
def get_all_tads(tadfile, chromosome):
    return get_tad(tadfile, chromosome, None)


#########################################
def get_alignment(species, assembly, chrom, region_from, region_to, scratch, outfile):
    # sudo apt install phast.v1_5.x86_64.deb (contains util maf_parse)
    maf_region_extraction_tool = "/usr/bin/maf_parse"
    # mafs come from here http://hgdownload.cse.ucsc.edu/downloads.html
    # http://hgdownload.cse.ucsc.edu/goldenPath/hg19/multiz100way/
    # or whichever species or assembly appropriate
    maf_file = "/storage/databases/ucsc/mafs/{}/{}/chr{}.maf".format(species, assembly, chrom)
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

    tmp_out = "{}/{}.maf".format(scratch, os.getpid())
    cmd = "{} -s {} -e {}  {} > {}".format(maf_region_extraction_tool,
                                           region_from, region_to, maf_file, tmp_out)
    subprocess.call(cmd, shell=True)

    cmd = "{} < {} > {}".format(maf2afa_tool, tmp_out, outfile)
    subprocess.call(cmd, shell=True)

    os.remove(tmp_out)


#########################################
if __name__ == "__main__":
    # print(ucsc_fragment_sequence('mm10',8, 57805369, 57805386))
    # get_alignment('mouse', 'mm10', 8, 57805369, 57805386, '/home/ivana/scratch', '/home/ivana/scratch/test.afa')
    print(ucsc_gene_coords('Hand2', "/storage/databases/ucsc/gene_ranges/mouse/mm9"))
