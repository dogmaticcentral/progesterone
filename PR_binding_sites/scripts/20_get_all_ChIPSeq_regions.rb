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

$random_check = false

##################################
def make_rand_file filename
     inf  = File.open "../data_raw/GSM857545_1_PR_oil_s_4_aligned.csv", "r"
     outf = File.open filename, "w"
     inf.each_line do |line|
          if line =~ /Sample/
               outf.write line
               next
          end
          sample, id, chrom, from, to = line.split "\t"
          from = from.to_i
          to = to.to_i
          if from>1.0e+6
               from -= 1.0e+6
               to -= 1.0e+6
          else
               from += 1.0e+6
               to += 1.0e+6
          end
          outf.write ["random", id, chrom, from, to].join("\t")
          outf.write("\n")
     end
     outf.close
     inf.close
end

##################################
region_sets = {}
if $random_check
     # we'll first make our random  intervals if they do not exist
     filename = "../data_raw/random.csv"
     make_rand_file filename if not File.exist? filename or File.zero? filename
     $random_chrom =  parse_chipseq_table filename
     region_sets   =  {"random"=>$random_chrom}
else
     # read in the regions detected in the presence of P4 and oil only
     $oil_chrom = parse_chipseq_table "../data_raw/GSM857545_1_PR_oil_s_4_aligned.csv"
     $p4_chrom  = parse_chipseq_table "../data_raw/GSM857546_2_PR_P4_s_1_aligned.csv"
     region_sets = {"oil"=> $oil_chrom, "p4"=> $p4_chrom}
     different_keys =  $p4_chrom.keys - $oil_chrom.keys
     if not different_keys.empty? then abort "oil and P4 do not seem to have the same chromosomes" end
end


region_sets.each do |label, region_set|
     region_set.sort_by {|k,v|  v.length}.each do |chrom, regions|
          file = File.open("../data_raw/#{label}_chr#{chrom}.fasta", "w")
          puts " chrom  #{chrom}  number of regions:   #{regions.length} "
          regions.each do |region|
               seq =  get_dna_region 'mm9', chrom, region.from, region.to
               file.write ">mm9_#{chrom}_#{region.from}_#{region.to}\n"
               (0...seq.length).step(50) { |i| file.write seq[i,50]+"\n"}
          end
          file.close
     end

end

