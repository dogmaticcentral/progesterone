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


##################################
# read in the regions detected in the presence of P4 and oil only
$oil_chrom = parse_chipseq_table "../data_raw/GSM857545_1_PR_oil_s_4_aligned.csv"
$p4_chrom  = parse_chipseq_table "../data_raw/GSM857546_2_PR_P4_s_1_aligned.csv"

different_keys =  $p4_chrom.keys - $oil_chrom.keys
if not different_keys.empty? then abort "oil and P4 do not seem to have the same chromosomes" end


$oil_chrom.sort_by {|k,v|  v.length}.each do |chrom, regions|
     file = File.open("../data_raw/oil_chr#{chrom}.fasta", "w")
     puts " chrom  #{chrom}  number of regions:   #{regions.length} "
     regions.each do |region|
          seq =  get_dna_region 'mm9', chrom, region.from, region.to
          file.write ">mm9_#{chrom}_#{region.from}_#{region.to}\n"
          (0..seq.length).step(35) { |i| file.write seq[i,50]+"\n"}
     end
     file.close
 end

$p4_chrom.sort_by {|k,v|  v.length}.each do |chrom, regions|
     file = File.open("../data_raw/p4_chr#{chrom}.fasta", "w")
     puts " chrom  #{chrom}  number of regions:   #{regions.length} "
     regions.each do |region|
          seq =  get_dna_region 'mm9', chrom, region.from, region.to
          file.write ">mm9_#{chrom}_#{region.from}_#{region.to}\n"
          (0..seq.length).step(35) { |i| file.write seq[i,50]+"\n"}
     end
     file.close
 end


