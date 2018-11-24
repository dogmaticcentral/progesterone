#!/usr/bin/python3
# mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
# -A skips auto rehash
from utils.mysqldb import *

#########################################
def main():
    # note the skip-auto-rehash option in .ucsc_myql_conf
    # it is the equivalent to -A on the mysql command line
    # means: no autocompletion, which makes mysql get up mych faster
    db     = connect_to_mysql("/home/ivana/.ucsc_mysql_conf")
    cursor = db.cursor()

    assembly = "hg18"
    #switch_to_db(cursor, "mm9") # mouse build name
    # no Y chromosome, we are looking at uterus tissue
    #chromosomes = ["chr"+str(x) for x in range(1,20)] + ["chrX"]
    switch_to_db(cursor, assembly) # human build name
    chromosomes = ["chr"+str(x) for x in range(1,23)] + ["chrX"]
    for chrom in chromosomes:
        print("downloading data for", chrom)
        outf = open("/storage/databases/ucsc/gene_ranges/human/{}/{}.csv".format(assembly,chrom), "w")
        outf.write( "\t".join(["name", "name2", "strand", "txStart", "txEnd"]) )
        qry  = "select name,  name2, strand, txStart, txEnd "
        qry += "from refGene "
        qry += "where chrom='%s' " % chrom
        qry += "and name like 'NM_%'"   # refseq says: NM_	mRNA	Protein-coding transcripts (usually curated)
        rows = search_db(cursor,qry)
        for row in rows:
            outf.write("\t".join( [ str(r) for r in row])+"\n")
        outf.close()
    cursor.close()
    db.close()

    
    
    return True


#########################################
if __name__ == '__main__':
    main()


