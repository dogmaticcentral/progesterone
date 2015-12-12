#!/usr/bin/env ruby

require_relative 'ruby_modules/region'
require_relative 'ruby_modules/parsers'
require_relative 'ruby_modules/utils'
require_relative 'ruby_modules/mysqlutils'
require_relative 'ruby_modules/dasutils'
include Parsers, Utils, MysqlUtils, DasUtils

$verbose = true
$scratch_space = "/tmp"
$extension = 25000

$testing = false

require 'bio' if $testing

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

# within 1kbp upstream from the first exon
# qry_genes = ['Epha4', 'Slc19a3', 'Lhx4', 'D930020B18Rik', 'Cd24a', 'Slc16a6', 'Gaa', 'Pik3ip1', 'Doc2b', 'Dbil5', 'Tob1', 'Ramp2', 'Smarcd2', 'Rin3', 'Pigh', 'Zfp410', 'Npc2', 'Ippk', 'Actn2', 'Hist1h2ac', 'Fbp1', 'Zbed3', 'Zmym5', 'Ephx2', 'Slc7a7', 'Irf9', 'C1qtnf6', 'Gga1', 'Asic1', 'Dlg1', 'Hcls1', 'Rcan1', 'Mpv17l', 'Ddr1', 'Mea1', '1700001C19Rik', 'Smad7', 'Cfl1', 'Shoc2', 'Aven', 'Tmem210', 'Fgf7', 'Glt28d2', 'S100a14', 'Rhbdl2', 'Ube4b', 'Lao1', 'Uchl1', 'Rest', 'Phc1', 'C1ra', 'Spsb2', 'Magohb', 'Sbk3', 'Gpr4', 'Cblc', 'Cebpg', 'Fzd4', 'Ccl25', 'Dctn6', 'Dpep1', 'Zfp599', 'Hspb2', 'Stag2', 'Fam46d', 'Nxt2', 'Amer1', 'Morf4l2' ]

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

# report for each gene (this should go)
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
          
          brackets = uniq_for_p4.map {|r| [r.from- min_tx_start, r.to- min_tx_start]}
          puts place_pr_regions brackets, strand, exons

          puts
          puts "binding motifs?"
          # http://press.endocrine.org/doi/10.1210/mend.7.4.8388996?url_ver=\
          # Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%3dpubmed
          puts "looking for  RgNaca NRN tgtNcY where R = A|G, Y = T|C, N = any"
          #                 cagaaca gtt tgttctg  2c7a
          nt2num = {'a' => 1, 't' => -1, 'g' => 2, 'c' =>-2 }
          uniq_for_p4.each do |region|
               puts "PR binding region #{region.from - min_tx_start} to #{region.to- min_tx_start} "
               seq =  get_dna_region 'mm9', chrom, region.from, region.to
               numeric = seq.split(//).map {|i| nt2num[i]}

               puts "   palindromic binding sites:"
               score = {}
               (7...numeric.length-7).each do |index|
                    # this is not true in the only cystal structure that we have
                    #sum  = numeric[index]>0 ? 0 : 1
                    sum = 0
                    sum += numeric[index-7]>0 ? 0 : 1
                    sum += numeric[index+7]<0 ? 0 : 1
                    2.upto(7).each {|j| sum +=  (numeric[index+j] + numeric[index-j])**2  }
                    score [ [index,seq[index-7,15] ] ]  = sum
               end
               score.sort_by{|k,v| v}.each { |k,v| puts "\t #{k[0]} #{k[1][0,6]} #{k[1][6,3]} #{k[1][9,6]}     #{v}"  if v <5 }

               puts "   half sites "
               score = {}
               (7...numeric.length-7).each do |index|
                    
                    #sum  = numeric[index]>0 ? 0.5 : 0
                    sum  = 0
                    sum  += 1 if seq[index-7,6] =~ /[ag]g[atcg]aca/
                    sum  += 1 if seq[index+2,6] =~ /tgt[atcg]c[tc]/
                    score [ [index,seq[index-7,15] ] ]  = sum
               end
               score.sort_by{|k,v| -v}.each { |k,v| puts "\t #{k[0]} #{k[1][0,6]} #{k[1][6,3]} #{k[1][9,6]}     #{v}"  if v > 0 }
               
              
          end
      
     end


 end

# cagaaca gtt tgttctg
