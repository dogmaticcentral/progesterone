#!/usr/bin/env ruby

require_relative 'ruby_modules/region'
require_relative 'ruby_modules/parsers'
require_relative 'ruby_modules/utils'
require_relative 'ruby_modules/mysqlutils'
require_relative 'ruby_modules/httputils'
include Parsers, Utils, MysqlUtils, HttpUtils
include Math

class ::String
     def matches? other_string, pattern
          return false if self.length != other_string.length
          match = true
          this_arr = self.split(//)
          other_arr = other_string.split(//)
          # this is just so awful ....
          (0 .. self.length).each {|i| (match &= (this_arr[i] == other_arr[i]) )  if pattern[i]}
          return match
     end
end

class ::Hash
  def has_pattern_key? (qry, pattern)
       self.keys.each do |k|
            return k if  k.matches?(qry,pattern)
       end
       return nil
  end
end

$verbose = true
$scratch_space = "/tmp"

$testing = false

require 'bio' if $testing

##################################
# read in the regions detected in the presence of P4 and oil only
$oil_chrom = parse_chipseq_table "../data_raw/GSM857545_1_PR_oil_s_4_aligned.csv"
$p4_chrom  = parse_chipseq_table "../data_raw/GSM857546_2_PR_P4_s_1_aligned.csv"

different_keys =  $p4_chrom.keys - $oil_chrom.keys
if not different_keys.empty? then abort "oil and P4 do not seem to have the same chromosomes" end

winlen  = 6
pattern = [true, true, true, false, true, true]
pattern = [true, true, false, true, true, true]

raw_counts = {}
uniq_containers = {}
number_uniq_containers = {}
parent_key = {}

$oil_chrom.sort_by {|k,v|  v.length}.each do |chrom, regions|
     puts " chrom  #{chrom}  number of regions:   #{regions.length} "

     regions.each do |region|
          #next if region.length > 500
          if region.length > 500
              start_pos = region.length/4.0
              end_pos = 3*region.length/4.0
          else
              start_pos = 0
              end_pos   = region.length
          end
          #puts "\t  #{region.from}  #{region.to} "
          seq =  get_dna_region 'mm9', chrom, region.from, region.to
          (start_pos ... end_pos-winlen).each do |index|
               frag = seq[index,winlen]
               if not raw_counts.has_key? frag
                    raw_counts[frag] = 0
                    number_uniq_containers[frag] = 0
                    uniq_containers[frag] = []
                    
               end
               raw_counts[frag]  += 1
               if not uniq_containers[frag].include? region
                    uniq_containers[frag].push region
                    number_uniq_containers[frag] += 1
               end
               
          end
          
     end
 end

#number_uniq_containers.sort_by {|k,v|  v}.each do |frag, count|
#     puts " #{frag}  #{count}  #{raw_counts[frag]} " #if count > 100
#end


# regroup the motifs that differ in the last two positions and cover the largesr portion of fragments
super_counts = {}
super_uniq_containers = {}
super_number_uniq_containers = {}

raw_counts.keys.each do |frag|

 
     key_match =  super_counts.has_pattern_key? frag, pattern 
     if not key_match
          super_counts[frag] = 0
          super_uniq_containers[frag] = []
          super_number_uniq_containers[frag] = 0
          key_match = frag
     end
          
     parent_key[frag] = key_match
     uniq_containers[frag].each  { |region|  (super_uniq_containers[key_match].push region) if not (super_uniq_containers[key_match].include? region) }
     
end

super_uniq_containers.sort_by {|k, v|  v.length}.each do |key_match, uniq_regions|
     
     print  " #{key_match}    #{uniq_regions.length}    "
     parent_key.select {|k2,v2| v2 == key_match}.keys.each {|frag| print frag, "  ", raw_counts[frag], "      "}
     puts

end
