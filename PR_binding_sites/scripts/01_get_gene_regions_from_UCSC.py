#!/usr/bin/python

import MySQLdb

########
def connect_to_mysql (conf_file):
    try:
        mysql_conn_handle = MySQLdb.connect(read_default_file=conf_file)
    except  MySQLdb.Error, e:
        print "Error connecting to mysql (%s) " % (e.args[1])
        sys.exit(1) 
    return mysql_conn_handle

########
def switch_to_db(cursor, db_name):
    qry = "use %s" % db_name
    rows = search_db(cursor, qry, verbose=False)
    if (rows):
        print rows
        return False
    return True


#######
def search_db(cursor, qry, verbose=False):
    try:
        cursor.execute(qry)
    except MySQLdb.Error, e:
        if verbose:
            print "Error running cursor.execute() for  qry: %s: %s " % (qry, e.args[1])
        return ["ERROR: " + e.args[1]]

    try:
        rows = cursor.fetchall()
    except MySQLdb.Error, e:
        if verbose:
            print "Error running cursor.fetchall() for  qry: %s: %s " % (qry, e.args[1])
        return ["ERROR: " + e.args[1]]

    if (len(rows) == 0):
        if verbose:
            print "No return for query %s" % qry
        return False

    return rows



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
        outf = open("gene_ranges."+chrom+".csv", "w")
        print  >>outf,  "\t".join( ["name", "name2", "txStart", "txEnd"] )
        qry  = "select g.name,  g.name2, g.txStart, g.txEnd "
        qry += "from refGene as g, gbCdnaInfo as i, refSeqStatus as s "
        qry += "where chrom='%s' " % chrom
        qry += "and i.acc=g.name  and i.type='mRNA' and i.mol='mRNA' "
        qry += "and s.mrnaAcc=g.name and s.status = 'Validated' "
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


