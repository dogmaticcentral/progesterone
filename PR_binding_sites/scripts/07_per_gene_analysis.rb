#!/usr/bin/env ruby

require_relative 'ruby_modules/region'
require_relative 'ruby_modules/parsers'
require_relative 'ruby_modules/utils'
require_relative 'ruby_modules/mysqlutils'
include Parsers, Utils, MysqlUtils

$verbose = true
$scratch_space = "/tmp"
$extension = 25000

##################################
# read in the regions detected in the presence of P4 and oil only
$oil_chrom = parse_chipseq_table "../data_raw/GSM857545_1_PR_oil_s_4_aligned.csv"
$p4_chrom  = parse_chipseq_table "../data_raw/GSM857546_2_PR_P4_s_1_aligned.csv"

different_keys =  $p4_chrom.keys - $oil_chrom.keys
if not different_keys.empty? then abort "oil and P4 do not seem to have the same chromosomes" end

# find ovelaps per chromosome
# uniq_for_p4 = {}
# $p4_chrom.each do |chrom,regions|
#      overlaps, uniq_for_p4[chrom] = find_overlapping_regions regions, $oil_chrom[chrom]
# end

# connect to ucsc, mm9 database
connection_handle = connect_to_mysql('/Users/ivana/.ucsc_myql_conf')
connection_handle.select_db ('mm9')

# find genes of interest
qry_genes = ['Fgf9', 'Spp1', 'Fgf2',  'Ptgs2', 'Hand2']
gene_coordinates = {}
qry_genes.each do |gene|
     gene_coordinates[gene] = {}
     qry = "select name, name2, chrom, strand, txStart, txEnd, exonStarts, exonEnds from refGene where name2='#{gene}'"
     connection_handle.query(qry).each do |row|
          gene_coordinates[gene][row['name']] = {
               chrom: row['chrom'].sub('chr',''),
               strand: row['strand'],
               region: Region.new(row['txStart'].to_i, row['txEnd'].to_i),
               exonStarts: row['exonStarts'],
               exonEnds: row['exonEnds']
          }                
     end
end
connection_handle.close

# are there any PR binding regions in the neoghborhood
qry_genes.each do |gene|
     min_tx_start = gene_coordinates[gene].map { |splice, coords| coords[:region].from }.min
     max_tx_end = gene_coordinates[gene].map { |splice, coords| coords[:region].to }.max
     translation_span = Region.new(min_tx_start, max_tx_end)
     covering_region = Region.new(min_tx_start-$extension, max_tx_end+$extension)

     chroms =  gene_coordinates[gene].map { |splice, coords| coords[:chrom]}.uniq
     (chroms.length == 1)  or abort "different chromosomes for different splices (?) for gene #{gene}"
     chrom = chroms[0]
     
     strands =  gene_coordinates[gene].map { |splice, coords| coords[:strand]}.uniq
     (strands.length == 1)  or abort "different strands for different splices (?) for gene #{gene}"
     strand = strands[0]
     
     overlaping_p4  = $p4_chrom[chrom].select {|region| region==covering_region}
     overlaping_oil = $oil_chrom[chrom].select {|region| region==covering_region}

     puts
     puts "="*30
     f = translation_span.from
     t = translation_span.to
     puts "gene #{gene}  chrom #{chrom}  strand #{strand}  from #{f}  to #{t}  length #{t-f} bp"
     puts "[in  the following all positions are given as measured from the earliest translation start (#{min_tx_start})]"
     
     write_exons = false
     if overlaping_p4.length > 0 
          puts "PR binding regions in the presence of progesterone"
          overlaping_p4.each do |region| 
               printf " %6d  %6d \n",  region.from - min_tx_start, region.to - min_tx_start
               write_exons |=  (region < translation_span)
          end
     else
          puts "no PR binding regions in the presence of progesterone "
     end
     if overlaping_oil.length > 0 
          puts "PR binding regions in the presence of vehicle (oil)"
          overlaping_oil.each do  |region|
               printf " %6d  %6d \n",  region.from - min_tx_start, region.to - min_tx_start
               write_exons |=  (region < translation_span)
          end
         
     else
          puts "no PR binding regions in the presence of vehicle (oil)"
     end

     overlaping_p4.length > 0 or next;
             
     overlaps, uniq_for_p4 = find_overlapping_regions overlaping_p4, overlaping_oil

     if overlaps.length > 0
          puts "note the overlapping regions between PR binding sites in the presence of progesterone and vehicle only: "
          overlaps.each do |ovlp|
               p4_region, oil_regions = ovlp
               printf " p4:  %6d  %6d \n",  p4_region.from - min_tx_start, p4_region.to - min_tx_start
               oil_regions.each do |oil_region|
                    printf " p4:  %6d  %6d \n",  oil_region.from - min_tx_start, oil_region.to - min_tx_start
               end
          end
     end

     
     uniq_for_p4.length > 0 or next;
     
     gene_coordinates[gene].each do |splice, coords|
          exon_starts = coords[:exonStarts].split(',').map {|i| i.to_i-min_tx_start}
          exon_ends = coords[:exonEnds].split(',').map {|i| i.to_i-min_tx_start}
          if write_exons
         
               puts "\nsplice: #{splice}"
               puts "region: from #{coords[:region].from- min_tx_start} to #{coords[:region].to- min_tx_start} "
               print "exonStarts "
               exon_starts.each {|i| printf "%6d  ", i}
               puts
               print "exonEnds   "
               exon_ends.each {|i| printf "%6d  ", i}
               puts
          end
          exons = exon_starts.zip(exon_ends)                   
          regions = uniq_for_p4.map {|r| [r.from- min_tx_start, r.to- min_tx_start]}
          puts place_pr_regions regions, strand, exons
      
     end

     puts
     # getting the seqeunce from das server
     #http://genome.ucsc.edu/cgi-bin/das/hg19/dna?segment=chr1:100000,200000
     # check on some exon seqs that we are getting the right seq back

 end
