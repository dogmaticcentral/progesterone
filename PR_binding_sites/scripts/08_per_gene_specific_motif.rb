#!/usr/bin/env ruby

require_relative 'ruby_modules/region'
require_relative 'ruby_modules/parsers'
require_relative 'ruby_modules/utils'
require_relative 'ruby_modules/mysqlutils'
require_relative 'ruby_modules/httputils'
include Parsers, Utils, MysqlUtils, HttpUtils
include Math

$verbose = true
$scratch_space = "/tmp"
$extension = 50000
$motif = 'tgttc'

$testing = false

require 'bio' if $testing


##################################
def score_region (region, chrom, strand, min_tx_start, max_tx_end, comment)
     score = {}
     seq =  get_dna_region 'mm9', chrom, region.from, region.to
     (0...seq.length-$motif.length).each do |index|
          
          sum  = 0
          subseq = seq[index,$motif.length]
          next if subseq!=$motif

         
          dist =  strand=='+' ? min_tx_start - region.to :  region.from - max_tx_end
          #next if dist < 0 or dist > 25000
          next if dist > $extension
          
          sum += 5*exp(- dist.abs/5000.0) 
          addr_from = region.from + index 
          addr_to = region.from   + index + $motif.length - 1
          score [ [index, addr_from, addr_to, subseq, comment ] ] = sum
          
     end
     return score
end

##################################
# read in the regions detected in the presence of P4 and oil only
$oil_chrom = parse_chipseq_table "../data_raw/GSM857545_1_PR_oil_s_4_aligned.csv"
$p4_chrom  = parse_chipseq_table "../data_raw/GSM857546_2_PR_P4_s_1_aligned.csv"

different_keys =  $p4_chrom.keys - $oil_chrom.keys
if not different_keys.empty? then abort "oil and P4 do not seem to have the same chromosomes" end
    
##################################
# connect to ucsc, mm9 database
connection_handle = connect_to_mysql('/Users/ivana/.ucsc_mysql_conf')
connection_handle.select_db ('mm9')

# find genes of interest
qry_genes = ['Scgb1a1']
    
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

# report for each gene 
qry_genes.each do |gene|
    
     min_tx_start = gene_coordinates[gene].map { |splice, coords| coords[:region].from }.min
     max_tx_end = gene_coordinates[gene].map { |splice, coords| coords[:region].to }.max
     translation_span = Region.new(min_tx_start, max_tx_end)
     covering_region  = Region.new(min_tx_start-$extension, max_tx_end+$extension)
     
     chroms =  gene_coordinates[gene].map { |splice, coords| coords[:chrom]}.uniq
     (chroms.length == 1)  or abort "different chromosomes for different splices (?) for gene #{gene}"
     chrom = chroms[0]
     
     strands =  gene_coordinates[gene].map { |splice, coords| coords[:strand]}.uniq
     (strands.length == 1)  or abort "different strands for different splices (?) for gene #{gene}"
     strand = strands[0]

     # are there any PR binding regions in the neoghborhood
     overlaping_p4  = $p4_chrom[chrom].select {|region| region==covering_region}
     overlaping_oil = $oil_chrom[chrom].select {|region| region==covering_region}

     puts
     puts "="*30
     f = translation_span.from
     t = translation_span.to
     puts "gene #{gene}  chrom #{chrom}  strand #{strand}  from #{f}  to #{t}  length #{t-f} bp"
     puts "[in  the following all positions are given as measured from the earliest translation start (#{min_tx_start})]"

     # report
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

     uniq_for_oil = []
     uniq_for_p4  = []
     overlaps, uniq_for_oil = find_overlapping_regions overlaping_oil, overlaping_p4

     if overlaps.length > 0
          puts "note the overlapping regions between PR binding sites in the presence of progesterone and vehicle only: "
          overlaps.each do |ovlp|
               oil_region, p4_regions = ovlp
               printf " oil:  %6d  %6d \n",  oil_region.from - min_tx_start, oil_region.to - min_tx_start
               p4_regions.each do |p4_region|
                    printf "  p4:  %6d  %6d \n",  p4_region.from - min_tx_start, p4_region.to - min_tx_start
               end
          end
     else
          overlaps, uniq_for_p4 = find_overlapping_regions overlaping_p4, overlaping_oil
     end
     

     #uniq_for_p4.length > 0 or next;
     
     gene_coordinates[gene].each do |splice, coords|
          exon_starts = coords[:exonStarts].split(',').map {|i| i.to_i-min_tx_start}
          exon_ends = coords[:exonEnds].split(',').map {|i| i.to_i-min_tx_start}
          exons = exon_starts.zip(exon_ends)                   
         
          if write_exons
         
               puts "\nsplice: #{splice}"
               puts "region: from #{coords[:region].from- min_tx_start} to #{coords[:region].to- min_tx_start} "
               print "exonStarts "
               exon_starts.each {|i| printf "%6d  ", i}
               puts
               print "exonEnds   "
               exon_ends.each {|i| printf "%6d  ", i}
               puts
               if $testing
                    puts "getting the exon seqs ..."
                    exons.each do |exon|
                         seq =  get_dna_region 'mm9', chrom, exon[0]+min_tx_start, exon[1]+min_tx_start
                         s = Bio::Sequence::NA.new(seq)
                         puts
                         puts seq
                         puts s.translate
                         puts s.translate (2)
                         puts s.translate (3)
                    end
              end
          end
          
         puts
         puts "binding motifs?"
         # http://press.endocrine.org/doi/10.1210/mend.7.4.8388996?url_ver=\
         # Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%3dpubmed
         puts "looking for #{$motif}\n"
         #                 cagaaca gtt tgttctg  2c7a
         score = {}
         overlaps.each do |ovlp|
             region, p4_regions = ovlp
             score.merge!(score_region region,chrom,strand,min_tx_start,max_tx_end,"oil_and_p4")
         end
         uniq_for_oil.each do |region|
               score.merge!(score_region region,chrom,strand,min_tx_start,max_tx_end,"oil_only")
          end
          uniq_for_p4.each do |region|
               score.merge!(score_region region,chrom,strand,min_tx_start,max_tx_end,"p4_only")
          end
          score.sort_by{|k,v| -v}.each { |k,v| puts "\t #{format("%4d",k[0])}   #{format("%8d",k[1])}  #{format("%8d",k[2])}  #{k[3][0,6]} #{k[3][6,3]} #{k[3][9,6]}  #{format("%15s",k[4])}  #{format("%.2f",v)} " }
          
     end
end

