#!/usr/bin/env ruby

require_relative 'ruby_modules/region'
require_relative 'ruby_modules/parsers'
require_relative 'ruby_modules/utils'
require_relative 'ruby_modules/httputils'
include Parsers, Utils,  HttpUtils
include Math

$verbose = false
$scratch_space = "/tmp"
$extension = 15000

$testing = false

require 'bio' if $testing

##################################
# read in the regions detected in the presence of P4 and oil only
$oil_chrom = parse_chipseq_table "../data_raw/GSM857545_1_PR_oil_s_4_aligned.csv"
$p4_chrom  = parse_chipseq_table "../data_raw/GSM857546_2_PR_P4_s_1_aligned.csv"

different_keys =  $p4_chrom.keys - $oil_chrom.keys
if not different_keys.empty? then abort "oil and P4 do not seem to have the same chromosomes" end


# find genes of interest
chromosomes = (1..19).map { |i| i.to_s} + ['X']
chromosomes.reverse!

# find the seqs most simlar to the PR binding motif ( RgNaca NRN tgtNcY where R = A|G, Y = T|C, N = any)
score = {}
chromosomes.each do |chrom|

     oil_and_p4_overlaps, uniq_for_p4 = find_overlapping_regions  $p4_chrom[chrom], $oil_chrom[chrom]
     oil_and_p4_binding_regions  = oil_and_p4_overlaps.map { |ovlp| ovlp[0] }
     
     gene_data, gene_ranges  = parse_gene_table(table_name:"../data_raw/gene_ranges.chr"+chrom+ ".csv",  extension: $extension, merge_splices: true)
     overlaps, genes_not_hit = find_overlapping_regions gene_ranges, uniq_for_p4

     # each overlap is a pair gene_region, pr_ranges
     puts "chrom: #{chrom}   genes hit: #{overlaps.length}"

     overlaps.each do |gene_region, pr_regions|
          # the first element of each splice_data array is txStart
          origin       = gene_region.from + $extension
          rev_origin   = gene_region.to   - $extension
          data         = gene_data[gene_region]
          gene_name    = data[0]
          splice_names = data[1]
          strand       = data[2]
          
          if $verbose
               puts "="*20
               puts " #{gene_name}    starts at: #{origin}    ends at #{rev_origin}"
               puts " splices:  #{splice_names}"
          end
          pr_regions.each do |pr_region|
               dist =  strand=='+' ? origin - pr_region.to  :  pr_region.from - rev_origin
               next if  dist < 0  or dist > 15000
             
               max_key = nil
               max_score = -1
               seq =  get_dna_region 'mm9', chrom, pr_region.from, pr_region.to
               (7...seq.length-7).each do |index|
                    
                    sum  = 0
                    subseq = seq[index-7,15]
                    sum += 1 if subseq[0] =~ /[ag]/
                    sum += 2 if subseq[1] == 'g'
                    if   subseq[3,3] =~ /aca/
                         sum += 4
                    elsif subseq[3,3] =~ /[ag][ct][ag]/
                         sum += 2
                    end
                    if   subseq[9,3] =~ /tgt/
                         sum += 4
                    elsif subseq[9,3] =~ /[ct][ag][ct]/
                         sum += 2
                    end
                    sum += 2 if subseq[13] == 'c'
                    sum += 1 if subseq[14] =~ /[ct]/
                    
                    dist =  strand=='+' ? origin - pr_region.to :  pr_region.from - rev_origin
                    #sum += 5*exp(- dist/500.0) end
                    region_length = pr_region.to - pr_region.from
                    next if index < region_length/4 or index > 3*region_length/4
                    sum += 5*exp(-dist/5000)
                    if sum > max_score
                         max_score = sum
                         key_string  = sprintf "%15s  %1s %2s  %10d  %10d ", gene_name, strand, chrom, origin, rev_origin
                         key_string += sprintf "  %6d  %4d  %3d    ", dist, region_length, index
                         key_string += sprintf "  %6s %3s %6s ", subseq[0,6], subseq[6,3],  subseq[9,6]
                         max_key = key_string
                    end

               end
               score [ max_key]  = max_score
          end
     end
     score.sort_by{|k,v| -v}.each { |k,v| puts "\t #{k}   #{format("%.2f",v)}"  if v > 5}
 end

outf = open("top_genes_half_site", 'w')
score.sort_by{|k,v| -v}.each { |k,v| outf.write("\t  #{k}   #{format("%.2f",v)}\n")  if v > 5 }
outf.close
 
# cagaaca gtt tgttctg
