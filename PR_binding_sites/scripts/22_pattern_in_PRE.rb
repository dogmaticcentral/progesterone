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

#######################################################################################
#######################################################################################
# define all those nice little aux fns we'll need below

@all_nts  = ['a', 'g', 't', 'c']
@purines  = ['a', 'g']
@pyrimidines  = ['a', 'g']

def build_motifs  alphabet, pattern, motif, ret_storage
     if pattern.length == 0
          ret_storage.push motif
     else
          if pattern[0]== :exact # exact match
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
 


def find_seq  dna_seq,  region_sets, region     
     region_sets.each do |label, region_set|
          region_set.sort_by {|k,v|  v.length}.each do |chrom, regions|
               if regions.include? region 
                    return dna_seq[label][chrom][ "mm9_#{chrom}_#{region.from}_#{region.to}"]
               end
          end
     end
     return nil
end

#######################################################################################
#######################################################################################
# start here, by reading in the regions of interest

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

     region_sets = {"p4"=> $p4_chrom}
     #region_sets = {"oil"=> $oil_chrom}
end


#######################################################################################
pattern = [:exact, :exact, :exact, :any, :exact,  :same_class_as_prev]
winlen  = pattern.length

motifs = []
build_motifs(@all_nts, pattern,'', motifs)
#
regex = {} # this is slow as hell, looks like I should have done it in perl
# to each motif assign  one regex for itself and pone for the complement - we'll consider all of these as equally valid
motifs.each do |motif|
     rev_compl = reverse_complement(motif)
     m1 = motif.sub("r","[ag]").sub("y", "[ct]")
     m2 = rev_compl.sub("r","[ag]").sub("y", "[ct]")
     regex[motif] =  Regexp.new(m1+"|"+m2 )
end

# since the regex search is so slow, now do something terrible: list all motifs that match the regex
seq2mot = {}
seq2motif @all_nts, pattern.length, '', regex , seq2mot

                                             

chromosomes = (1..19).map { |i| i.to_s} + ['X']
dna_seq = {}
if $random_check
     labels = ['random']
else
     labels = ['oil','p4']
end
labels.each do |label|
     dna_seq[label] = {}
     chromosomes.each do |chrom|
          dna_seq[label][chrom] = read_dna label, chrom
     end
end


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
          regions.each do |region|
               
               #next if region.length > 500
               total_regions += 1
               start_pos = (region.length/5.0).to_i
               end_pos   = (4*region.length/5.0).to_i

               seq = dna_seq[label][chrom][ "mm9_#{chrom}_#{region.from}_#{region.to}"]
               next if seq =~ /n/ # I have 4 cases of this altogether
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
     next if ($random_check and pct < 70) or (!$random_check and pct < 80)
     #next if  motif != "tgg.cy"
     print  " #{motif}    #{uniq_regions.length}/#{total_regions}  (#{format("%.0f",pct)}%)    "
     # how many of these have the reverse complement back to back?
     m1 = motif.sub("r","[ag]").sub("y", "[ct]")
     m2 = reverse_complement(motif).sub("r","[ag]").sub("y", "[ct]")
     pattern = m2+"..."+m1
     #pattern = m1+"..."+m2
     reggie =  Regexp.new(pattern)
     hits = 0
     uniq_regions.each do |region|
          seq = find_seq  dna_seq, region_sets, region
          next if not seq
          next if seq =~ /n/
          match_data =  seq.scan(reggie)
          if match_data.length > 0             
               hits += 1
               #match_data.each do |matching_seq|
               #     puts "\t\t " + matching_seq
               #end
               #puts  "\t\t " + "-"*15
         end
     end
     pct = hits*100.0/uniq_regions.length
     puts "    #{pattern}   #{hits}/#{uniq_regions.length}  (#{format("%.1f",pct)}%) "
end









###############################################################
###############################################################
if $sanity_check
     # check
     motif = "tgt.cy"
     tot_found = 0
     tot_regions = 0
     region_sets.each do |label, region_set|
          region_set.sort_by {|k,v|  v.length}.each do |chrom, regions|
               puts " #{label} chrom  #{chrom}  number of regions:   #{regions.length} "
               regions.each do |region|
                    
                    #next if region.length > 1000
                    tot_regions += 1
                    start_pos = (region.length/5.0).to_i
                    end_pos   = (4*region.length/5.0).to_i

                    seq = dna_seq[label][chrom][ "mm9_#{chrom}_#{region.from}_#{region.to}"]
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
