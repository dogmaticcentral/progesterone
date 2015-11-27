#!/usr/bin/env ruby

class Region
     attr_accessor :from, :to
     def initialize from, to
          @from = from
          @to = to
     end

     def  overlaps? reg2
          not (@to < reg2.from or reg2.to < @from)
     end

     alias == overlaps?
end

# & is not applicable to this case, bcs two regions are equivalent if they
# overlap, not if they have the same from and to (so which one of them should be
# used as the representative in the intersection?)
def find_overlapping_regions regions_list_1, regions_list_2
     overlaps = Array.new
     uniq_for_list_1 = Array.new
     regions_list_1.each do |region_1|
          overlap = regions_list_2.select { |region_2|  region_2  == region_1 }
          overlap.empty? ? uniq_for_list_1.push(overlap) :  overlaps.push( [region_1, overlap])
     end
     return overlaps, uniq_for_list_1
end

def parse_table table_name
     chromosome = Hash.new
     File.readlines(table_name).drop(1).each do |line|
          sample, id, chrom, from, to = line.split "\t"
          (chromosome.key? chrom)  || chromosome[chrom] = Array.new 
          chromosome[chrom].push(Region.new(from.to_i,to.to_i))
     end
     return chromosome
end

$oil_chrom = parse_table "GSM857545_1_PR_oil_s_4_aligned.csv"
$p4_chrom  = parse_table "GSM857546_2_PR_P4_s_1_aligned.csv"

different_keys =  $p4_chrom.keys - $oil_chrom.keys
if not different_keys.empty? then abort "oil and P4 do not seem to have the same chromosomes" end


def find_overlaps_per_chromosome
     total = 0
     $p4_chrom.each do |chrom,regions|
          overlaps, uniq_for_p4 = find_overlapping_regions regions, $oil_chrom[chrom]
          total += uniq_for_p4.length
          puts "chrom #{chrom}  number of overlapping regions #{overlaps.length} "
          puts "                uniq to P4 #{uniq_for_p4.length}  "
     end
     puts total
end


find_overlaps_per_chromosome

