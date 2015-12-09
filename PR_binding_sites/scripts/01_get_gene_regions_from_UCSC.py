#!/usr/bin/python

from   python_modules.mysqldb import *

#########################################
def main():
    # note the skip-auto-rehash option in .ucsc_myql_conf
    # it is the equivalent to -A on the mysql command line
    # means: no autocompletion, which makes mysql get up mych faster
    db     = connect_to_mysql("/Users/ivana/.ucsc_myql_conf")
    cursor = db.cursor()
    
    switch_to_db(cursor, "mm9") # mouse build name
    # no Y chromosome, we are looking at uterus tissue
    chromosomes = ["chr"+str(x) for x in range(1,20)] + ["chrX"]
    for chrom in chromosomes:
        print "downloading data for", chrom
        outf = open("../data_raw/gene_ranges."+chrom+".csv", "w")
        print  >>outf,  "\t".join( ["name", "name2", "txStart", "txEnd"] )
        qry  = "select g.name,  g.name2, g.txStart, g.txEnd "
        qry += "from refGene as g, gbCdnaInfo as i, refSeqStatus as s "
        qry += "where chrom='%s' " % chrom
        qry += "and i.acc=g.name  and i.type='mRNA' and i.mol='mRNA' "
        qry += "and s.mrnaAcc=g.name and s.status in ('Validated', 'Reviewed') "
        rows = search_db(cursor,qry)
        
        for row in rows:
            print  >>outf,  "\t".join( [ str(r) for r in row] )
        outf.close()
    cursor.close()
    db.close()

    
    
    return True


#########################################
if __name__ == '__main__':
    main()


