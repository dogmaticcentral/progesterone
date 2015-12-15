#!/usr/bin/env ruby

require_relative 'ruby_modules/region'
require_relative 'ruby_modules/parsers'
require_relative 'ruby_modules/utils'
include Parsers, Utils

$verbose = false
$scratch_space = "/tmp"


##################################
# read in the regions detected in the presence of P4 and oil only
$oil_chrom = parse_chipseq_table "../data_raw/GSM857545_1_PR_oil_s_4_aligned.csv"
$p4_chrom  = parse_chipseq_table "../data_raw/GSM857546_2_PR_P4_s_1_aligned.csv"

different_keys =  $p4_chrom.keys - $oil_chrom.keys
if not different_keys.empty? then abort "oil and P4 do not seem to have the same chromosomes" end

# find ovelaps per chromosome
uniq_for_p4 = {}
$p4_chrom.each do |chrom,regions|
     overlaps, uniq_for_p4[chrom] = find_overlapping_regions regions, $oil_chrom[chrom]
end

# find overlaps between regions that pop up only in the presence of P4
chromosomes = (1..19).map { |i| i.to_s} + ['X']
total_genes_hit = 0
chromosomes.each do |chrom|
     
     gene_data, gene_ranges  = parse_gene_table "../data_raw/gene_ranges.chr"+chrom+ ".csv"
     overlaps, genes_not_hit = find_overlapping_regions gene_ranges, uniq_for_p4[chrom]

     genes_hit = {}
     
     overlaps.each do |overlap|
          gene_region, pr_ranges = overlap
          name, name2, strand, txStart, txEnd, exonStarts, exonEnds  = gene_data[gene_region].split "\t"
          if not genes_hit.key? name2
               genes_hit[name2] = { splices: [[name, strand, txStart.to_i, txEnd.to_i]],
                    pr_binding_sites:  pr_ranges}
          else
               genes_hit[name2][:splices].push [name, strand, txStart.to_i, txEnd.to_i]
               pr_ranges.each do |this_region|
                    ranges_seen = genes_hit[name2][:pr_binding_sites].select { |prb| prb === this_region }
                    if ranges_seen.empty? then genes_hit[name2][:pr_binding_sites].push(this_region) end
               end
          end
     end

     if $verbose
          puts "chrom #{chrom}  PR binding regions uniq to P4 #{uniq_for_p4[chrom].length}  "
          puts "         genes:  #{gene_data.length}  "
          puts "         genes hit:  #{genes_hit.keys.length}  "
     else
          outstr = ""
     end
     total_genes_hit += genes_hit.length

     genes_hit.each do |name, hash|
          # the first element of each splice_data array is txStart
          strand = hash[:splices][0][1]
          origin = hash[:splices].map { |splice_data| splice_data[2] }.min
          rev_origin = hash[:splices].map { |splice_data| splice_data[3] }.max
          if $verbose
               puts "="*20
               puts " #{name}    starts at: #{origin} "
               puts " splices: "
               hash[:splices].each do |splice_data|
                    splice_name, txStart, txEnd = splice_data
                    puts " \t #{splice_name}    #{txStart-origin}   #{txEnd-origin}  "
               end
               puts " PR binding sites: "
               hash[:pr_binding_sites].each do |pr_region|
                    puts " \t #{pr_region.from-origin}    #{pr_region.to-origin}  "
               end
          else
               outstr +=   "#{name}\t#{strand}\t#{origin}\t#{rev_origin}\t"
               outstr +=   hash[:splices].map { |splice_data| splice_data[0] }.join","
               outstr +=  "\t"
               outstr +=   hash[:pr_binding_sites].map {|bs|  "#{bs.from-origin}..#{bs.to-origin}" }.join","
               outstr +=  "\n"
          end
     end
     $verbose || File.write($scratch_space+'/p4_chr'+chrom+'.txt', outstr)
end     
puts total_genes_hit
