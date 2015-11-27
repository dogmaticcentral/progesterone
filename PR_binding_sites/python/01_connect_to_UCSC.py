#!/usr/bin/env python

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
    
    for chrom in range(1,20):
        qry = "select refGene.name, refGene.chrom, refGene.name2  from refGene, gbCdnaInfo "
        qry += "where chrom='chr%s' " % chrom
        qry += "and gbCdnaInfo.acc=refGene.name  and gbCdnaInfo.mol='mRNA' "
        qry += "limit 10"
        rows = search_db(cursor,qry)
        for row in rows[:10]:
            print row
 
    cursor.close()
    db.close()

    
    
    return True


#########################################
if __name__ == '__main__':
    main()


