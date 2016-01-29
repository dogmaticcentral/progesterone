#!/usr/bin/env ruby

# just by looking at the regions that the ChIPSeq pulls down - can we tell the  motif?
# it looks like tgt.cy and be found ins some 86% of random cases (random regions)
# vs 91% of cases labeled oil
# vs 94% of cases labeled p4
# vs 93% of cases if both o&p

require_relative 'ruby_modules/region'
require_relative 'ruby_modules/parsers'
require_relative 'ruby_modules/utils'
require_relative 'ruby_modules/mysqlutils'
require_relative 'ruby_modules/httputils'

require 'benchmark'

include Parsers, Utils, MysqlUtils, HttpUtils
include Math


$verbose = true
$scratch_space = "/tmp"
$sanity_check = false
$random_check = false

##################################
region_sets = {}
if $random_check
     $random_chrom =  parse_chipseq_table "../data_raw/random.csv"
     region_sets   =  {"random"=>$random_chrom}
else
     # read in the regions detected in the presence of P4 and oil only
     $oil_chrom = parse_chipseq_table "../data_raw/GSM857545_1_PR_oil_s_4_aligned.csv"
     $p4_chrom  = parse_chipseq_table "../data_raw/GSM857546_2_PR_P4_s_1_aligned.csv"
     region_sets = {"oil"=> $oil_chrom, "p4"=> $p4_chrom}
     different_keys =  $p4_chrom.keys - $oil_chrom.keys
     if not different_keys.empty? then abort "oil and P4 do not seem to have the same chromosomes" end
end
region_sets = {"p4"=> $p4_chrom}


@all_nts  = ['a', 'g', 't', 'c']
@purines  = ['a', 'g']
@pyrimidines  = ['a', 'g']

def build_motifs  alphabet, pattern, motif, ret_storage
     if pattern.length == 0
          ret_storage.push motif
     else
          if pattern[0]== :same # exact match
               alphabet.each { |c| build_motifs(alphabet, pattern[1..-1], motif+c, ret_storage ) }
          elsif pattern[0]==:same_class_as_prev # class here refers to purine vs pyrimidines
               if @purines.include? motif[-1]
                    build_motifs(@purines, pattern[1..-1], motif+'r', ret_storage)
               else
                    build_motifs(@pyrimidines, pattern[1..-1], motif+'y', ret_storage)
               end
          else
               build_motifs(alphabet, pattern[1..-1], motif+'.', ret_storage)
          end
     end
end

def seq2motif alphabet, length, parent_seq, regex, seq2mot
     if length == 0
          seq2mot[parent_seq] = []
          regex.each do  |motif,reg| 
             seq2mot[parent_seq].push(motif) if parent_seq =~ reg
          end
     else
          alphabet.each { |c| seq2motif alphabet, length-1, parent_seq+c, regex, seq2mot }
     end
end

def initialize_hash  alphabet, length, parent_seq, hash
     if length == 0
          hash[parent_seq] = 0
     else
          alphabet.each { |c| initialize_hash alphabet, length-1, parent_seq+c, hash }
     end
end
 

def reverse_complement blah
     compl = {'a'=>'t', 't'=>'a','c'=>'g',  'g'=>'c', '.'=>'.', 'r'=>'y',  'y'=>'r'}
     retstr = ''
     blah.downcase.reverse.split(//).each {|c|  retstr += compl[c]}
     return retstr
end

pattern = [:same, :same, :same, :any, :same,  :same_class_as_prev]
#pattern = [true, true, true, false, true]
winlen  = pattern.length

motifs = []
build_motifs(@all_nts, pattern,'', motifs)
#
regex = {} # this is slow as hell, looks like I should have done it in perl
# to each motif assig one regex for itself and pone for the complement - we'll consider all of these as equally valid
motifs.each do |motif|
     rev_compl = reverse_complement(motif)
     m1 = motif.sub("r","[ag]").sub("y", "[ct]")
     m2 = rev_compl.sub("r","[ag]").sub("y", "[ct]")
     regex[motif] =  Regexp.new(m1+"|"+m2 )
end

# since the regex search is so slow, now do something terrible: list all motifs that match the regex
seq2mot = {}
seq2motif @all_nts, pattern.length, '', regex , seq2mot

                                             
raw_counts = {}
uniq_containers = {}
motifs.each do |motif| 
     uniq_containers[motif] = []
end

raw_counts = {}
initialize_hash  @all_nts, winlen, '', raw_counts

total_regions = 0

region_sets.each do |label, region_set|
     region_set.sort_by {|k,v|  v.length}.each do |chrom, regions|
          puts " #{label} chrom  #{chrom}  number of regions:   #{regions.length} "
          dna_seq = read_dna label, chrom
          regions.each do |region|
               
               #next if region.length > 1000
               total_regions += 1
               start_pos = (region.length/5.0).to_i
               end_pos   = (4*region.length/5.0).to_i

               seq = dna_seq[ "mm9_#{chrom}_#{region.from}_#{region.to}"]
               next if seq =~ /n/
               motifs_in_this_region = []
               
               (start_pos ... end_pos-winlen).each do |index|
                    frag = seq[index,winlen]
                    raw_counts[frag]  += 1
                    motifs = seq2mot[ frag ]
                    # which motifs does this frag make happy
                    motifs.each do |motif| 
                         motifs_in_this_region.push(motif) if not  motifs_in_this_region.include?  motif
                    end
                    
               end

               motifs_in_this_region.each {|m| uniq_containers[m].push region}
               
          end
     end
end

mot2seq = {}
uniq_containers.each {|mot,conts| mot2seq[mot] = []}
seq2mot.each do |seq, motifs|
     motifs.each {|m| mot2seq[m].push(seq) }
end


uniq_containers.sort_by {|k, v|  v.length}.each do |motif, uniq_regions|
     pct = uniq_regions.length*100.0/total_regions
     next if pct < 80
     print  " #{motif}    #{uniq_regions.length}/#{total_regions}  (#{format("%.0f",pct)}%)    "
     #mot2seq[motif].each {|seq| print "  #{seq}  #{raw_counts[seq]} "}
     puts
          
end

###############################################################
if $sanity_check
     # check
     motif = "tgt.cy"
     tot_found = 0
     tot_regions = 0
     region_sets.each do |label, region_set|
          region_set.sort_by {|k,v|  v.length}.each do |chrom, regions|
               puts " #{label} chrom  #{chrom}  number of regions:   #{regions.length} "
               dna_seq = read_dna label, chrom
               regions.each do |region|
                    
                    #next if region.length > 1000
                    tot_regions += 1
                    start_pos = (region.length/5.0).to_i
                    end_pos   = (4*region.length/5.0).to_i

                    seq = dna_seq[ "mm9_#{chrom}_#{region.from}_#{region.to}"]
                    next if seq =~ /n/
                    found = false
                    (start_pos ... end_pos-winlen). each do |index|
                         frag = seq[index,winlen]
                         raw_counts[frag]  += 1
                         motifs = seq2mot[ frag ]
                         if motifs.include? motif
                              puts " #{index}/#{region.length}   #{frag} "
                              tot_found += 1 if  not found
                              found = true
                              
                         end
                    end

               end
          end
     end

     pct = tot_found*100.0/tot_regions
     puts  " #{motif}    #{tot_found}/#{tot_regions}  (#{format("%.0f",pct)}%) "
end
