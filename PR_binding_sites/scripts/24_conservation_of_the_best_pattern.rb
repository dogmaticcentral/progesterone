#!/usr/bin/env ruby

# just by looking at the regions that the ChIPSeq pulls down - can we tell the  motif?
# it looks like tgt.cy and be found ins some 86% of random cases (random regions)
# vs 91% of cases labeled oil
# vs 94% of cases labeled p4
# vs 93% of cases if both o&p

require_relative 'ruby_modules/region'
require_relative 'ruby_modules/parsers'
require_relative 'ruby_modules/utils'

include Parsers, Utils

require_relative 'ruby_modules/httputils'
include  HttpUtils


motif = "tgt.cy"

m1 = motif.sub("r","[ag]").sub("y", "[ct]")
m2 = reverse_complement(motif).sub("r","[ag]").sub("y", "[ct]")

pattern = m2+"..."+m1
regexp_double_site =  Regexp.new(pattern)



$oil_chrom = parse_chipseq_table "../data_raw/GSM857545_1_PR_oil_s_4_aligned.csv"
$p4_chrom  = parse_chipseq_table "../data_raw/GSM857546_2_PR_P4_s_1_aligned.csv"
region_sets = {"oil"=> $oil_chrom, "p4"=> $p4_chrom}
different_keys =  $p4_chrom.keys - $oil_chrom.keys
if not different_keys.empty? then abort "oil and P4 do not seem to have the same chromosomes" end

region_sets = {"p4"=> $p4_chrom}


labels = ['oil','p4']
chromosomes = (1..19).map { |i| i.to_s} + ['X']

dna_seq = {}
labels.each do |label|
     dna_seq[label] = {}
     chromosomes.each do |chrom|
          dna_seq[label][chrom] = read_dna label, chrom
     end
end


region_sets.each do |label, region_set|
     region_set.sort_by {|k,v|  v.length}.each do |chrom, regions|
          puts " #{label} chrom  #{chrom}  number of regions:   #{regions.length} "
          regions.each do |region|
               
               next if region.length > 500
               start_pos = (region.length/5.0).to_i
               end_pos   = (4*region.length/5.0).to_i

               seq = dna_seq[label][chrom][ "mm9_#{chrom}_#{region.from}_#{region.to}"]

               puts seq[0,100]
               seq2 =  get_dna_region 'mm9', chrom, region.from, region.to
               puts seq2[0,100]

               
               # I have no friggin idea how this works but it does
               positions = seq.enum_for(:scan, regexp_double_site).map { Regexp.last_match.begin(0) }.sort
               if positions.length > 0
                    puts  "\t  #{region.from}   #{region.to} \n"
                    positions.each do |i|
                         puts "\t\t   #{i}   #{seq[i,15]} "
                         puts "\t\t         #{seq2[i,15]} "
                         exit
                         alignment_file = "alignments/#{chrom}_#{region.from+i}_#{region.from+i+14}.afa"
                         if not File.exist? alignment_file
                              get_alignment alignment_file, chrom, region.from+i, region.from+i+14
                         end
                    end
                    exit
               end
               
           end
     end
end
