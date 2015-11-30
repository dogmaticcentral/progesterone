#!/usr/bin/python -u

from   python_modules.mysqldb import *

#########################################
def place_pr_regions (pr_regions, strand, exons):
    descr_str = ""
    regions = [ map ( lambda x: int(x), pair.split("..") ) for pair in pr_regions.split (",")]
    orig_regions = regions[:]
    if strand=='-':
        origin = exons[-1][1]
        regions = map (lambda x: [-(x[1]-origin), -(x[0]-origin)], regions)
        regions.reverse()
        exons = map (lambda x: [-(x[1]-origin), -(x[0]-origin)], exons)
        exons.reverse()
    
    first = True
    for region in regions:
        if first:
            first = False
        else:
            descr_str += "\t"
        if strand=='-':
            orig_index = len(regions) - regions.index(region) - 1
            orig_region = orig_regions[orig_index]
            descr_str += "%d..%d:   " % (orig_region[0], orig_region[1])
        else:
            descr_str += "%d..%d:   " % (region[0], region[1])
        if  region[1] < exons[0][0]:
            diff = exons[0][0]-region[1]
            if diff < 1000:
                descr_str += "%d bp upstream from the first exon" % diff
            else:
                descr_str += "%.1f kbp upstream from the first exon" % (diff/1000.0 )               
        elif  exons[-1][1] < region[0]:
            diff = region[0] -  exons[-1][1]
            if diff < 1000:
                descr_str += "%d bp downstream from the last exon" % diff
            else:
                descr_str += "%.1f kbp downstream from the last exon" % (diff/1000.0) 
        else:
            start = []
            if  region[0] >= exons[0][0]:
                for i in range(len(exons)):
                    if region[0] <= exons[i][1]:
                        start = ["exon", i+1]
                        break
                    if i+1 < len(exons) and region[0] < exons[i+1][0]:
                        start = ["intron", i+1]
                        break
            end = []
            if region[0] <= exons[-1][1]:
                for i in range(len(exons)):
                    if region[1] <= exons[i][1]:
                        end = ["exon", i+1]
                        break
                    if i+1 < len(exons) and region[1] < exons[i+1][0]:
                        end = ["intron", i+1]
                        break
            if len(start) ==0 and len(end)== 0:
                descr_str += "err locating pr region (?)"
            elif len(end)== 0:
                descr_str +=  "starts in  %s %d, ends after the last exon " % (start[0], start[1])
            elif len(start)== 0:
                descr_str +=  "starts before the first exon, ends in  %s %d " % (end[0], end[1])
            elif start[0] == end[0] and  start[1] == end[1]:
                descr_str += "within  %s %d " % (start[0], start[1])
            else:
                descr_str += "starts in  %s %d, ends in   %s %d" % (start[0], start[1], end[0], end[1])
                
    return descr_str
    

#########################################
def main():

    db     = connect_to_mysql("/Users/ivana/.ucsc_myql_conf")
    cursor = db.cursor()
    switch_to_db(cursor, "mm9") # mouse build name

    chromosomes = ["chr"+str(x) for x in range(1,20)] + ["chrX"]
    for chrom in chromosomes:
        infile  = open ('/tmp/p4_'+chrom+'.txt', "r")
        outfile = open ("../data_processed/annotated_p4_pr_regions_"+chrom+'.csv', "w")
        print >>outfile, "\t".join (["name", "splice",  "description", "strand",  "trsl start", 
                                     "trsl end", "length", "number of exons", "PR binding region (from trsl start)"]) 
        for line in infile:
            line = line.rstrip()
            [name, origin, rev_origin, splices_str, pr_regions] = line.split("\t")
            origin = int(origin)
            rev_origin = int (rev_origin)
            splices = splices_str.split (",")
            
            prev_id = ""
            exons  = {}
            description_ids = []
            for splice in splices:
                qry  = "select name2, strand, exonStarts, exonEnds "
                qry += "from  refGene  "
                qry += "where chrom='%s' " % chrom
                qry += " and name='%s' " % splice
                rows = search_db(cursor,qry)
                if len(rows) > 1:
                    print "splice name %s not uniq:  we shouldn't be here ..." % splice
                    next
                name2, strand,  exonStarts, exonEnds = rows[0]
                if name2 != name:
                    print "gene name mismatch: %s  %s: we shouldn't be here ..." % (name,  name2)
                    next

                if exonStarts[-1] == ",": exonStarts= exonStarts[:-1]
                if exonEnds[-1]   == ",": exonEnds= exonEnds[:-1]
                eS = [ int(x)-int(origin) for x in exonStarts.split (",")]
                eE = [ int(x)-int(origin) for x in exonEnds.split (",")]
                exons  =  zip (eS, eE)

                qry  = "select description.name from description, gbCdnaInfo "
                qry += "where gbCdnaInfo.acc ='%s' " % splice
                qry += "and gbCdnaInfo.description= description.id"
                rows = search_db(cursor,qry)
                if len(rows) > 1:
                    print "multiple descriptions: "
                    for row in rows: print row[0]
                    print "we shouldn't be here ..."
                    next
                description = rows[0][0]
                if 'Mus musculus' in description[:12]:
                    description = description[13:]
                pr_region_description_str = place_pr_regions (pr_regions, strand, exons)
                if strand=='-':
                    strand_string = "negative"
                else:
                    strand_string = "positive"
                print  >>outfile, "\t".join (map (lambda x: str(x), [name, splice,  description, strand_string,  origin, 
                                  rev_origin, rev_origin-origin, len(exons), pr_region_description_str]) )
        outfile.close()
        infile.close()
 
    cursor.close()
    db.close()

#########################################
if __name__ == '__main__':
    main()
