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



def build_motifs  alphabet, pattern, motif, ret_storage
     if pattern.length == 0
          ret_storage.push motif
     else
          alphabet.each { |c| build_motifs(alphabet, pattern[1..-1], motif+(pattern[0]?c:'.'), ret_storage ) }
         
     end
end

def reverse_complement blah
     compl = {'a'=>'t', 't'=>'a','c'=>'g',  'g'=>'c', '.'=>'.'}
     retstr = ''
     blah.downcase.reverse.split(//).each {|c|  retstr += compl[c]}
     return retstr
end

winlen  = 6
pattern = [true, true, true, false, true, true]
alphabet = ['a', 'g', 't', 'c']

motifs = []
build_motifs(alphabet,pattern,'', motifs)
regex = {}
# to each motif assig one regex for itself and pone for the complement - we'll consider all of these as equally valid
motifs.each { |motif|  regex[motif] =  Regexp.new(motif+"|"+reverse_complement(motif) )}
                                             
raw_counts = {}
uniq_containers = {}
motifs.each do |motif| 
     raw_counts[motif] = 0
     uniq_containers[motif] = []
end



total_regions = 0
$oil_chrom.sort_by {|k,v|  v.length}.each do |chrom, regions|
     puts " chrom  #{chrom}  number of regions:   #{regions.length} "
     regions.each do |region|
          #next if region.length > 500
          total_regions += 1
          #if region.length > 500
          #    start_pos = region.length/4.0
          #    end_pos = 3*region.length/4.0
          #else
          start_pos = 0
          end_pos   = region.length

          seq =  get_dna_region 'mm9', chrom, region.from, region.to
          (start_pos ... end_pos-winlen).each do |index|
               frag = seq[index,winlen]
               # which motifs does this frag make happy
               motifs.each do |motif| 
                    next if not regex[motif].match(frag)
                    raw_counts[motif]  += 1
                    if not uniq_containers[motif].include? region
                         uniq_containers[motif].push region
                    end
                end
              
          end
          
     end
 end



uniq_containers.sort_by {|k, v|  v.length}.each do |motif, uniq_regions|
     pct = uniq_regions.length*100.0/total_regions
     next if pct < 70
     puts  " #{motif}    #{uniq_regions.length}/#{total_regions}  (#{format("%.0f",pct)}%)    "

end
