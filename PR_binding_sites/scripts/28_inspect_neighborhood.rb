#!/usr/bin/env ruby  -W0

require_relative 'ruby_modules/parsers'
require_relative 'ruby_modules/mysqlutils'
require_relative 'ruby_modules/httputils'
require_relative 'ruby_modules/utils'

include Parsers, Utils, MysqlUtils, HttpUtils


chromosomes = (1..19).map { |i| i.to_s} + ['X']

motif_positions = {}
chromosomes.each {|chrom| motif_positions[chrom] = []}
     
File.readlines('best_alignments').each do |line|
     
     name, score, no_species = line.split "\s"
     chrom, start, stop = name.gsub("alignments/", "").gsub(".afa","").split "_"
     motif_positions[chrom].push ([start.to_i,stop.to_i])
end


motif = "tgt.cy"

m1 = motif.sub("r","[ag]").sub("y", "[ct]")
m2 = reverse_complement(motif).sub("r","[ag]").sub("y", "[ct]")

half_site_left  =  Regexp.new(m2)
half_site_right =  Regexp.new(m1)



chromosomes.each do |chrom|
     next if motif_positions[chrom].size == 0
     gene_data, gene_ranges  = parse_gene_table(table_name:"../data_raw/gene_ranges.chr"+chrom+ ".csv",extension:0, merge_splices:true)

     motif_positions[chrom].each do |m|
          motif_from, motif_to = m
          dist = {}
          strand = {}
          gene_data.each_value do |v|
               name, splice_names, str, tx_from, tx_to = v
               strand[name] = str
               if motif_to < tx_from
                    dist[name] = motif_to - tx_from 
               elsif tx_to < motif_from
                    dist[name] = motif_from - tx_to
               else
                    dist[name] = 0
               end
          end
          dist.sort_by {|k,v| v.abs}.to_h.keys[0,3].each do |name|
               next if dist[name].abs>20000
               next if dist[name] == 0
               next if strand[name]=="+" and  dist[name] >0
               next if strand[name]=="-" and  dist[name] <0
               puts "#{chrom}   #{motif_from} #{motif_to}   #{name}  #{strand[name]}    #{dist[name]}"
               # get the extended sequence
               nbrhood = 250
               seq =  get_dna_region 'mm9', chrom, motif_from-nbrhood, motif_to+nbrhood
               before = true
               [ seq[0...nbrhood], seq[-nbrhood..-1] ].each do |nbrhood_seq|
                    if before
                         puts "nbrhood left"
                         before = false
                    else
                         puts "nbrhood right"
                    end
                    [ half_site_left, half_site_right].each do |half_site|
                         positions = nbrhood_seq.enum_for(:scan, half_site).map { Regexp.last_match.begin(0) }.sort
                         puts "\thalf sites: #{positions.length}"
                         positions.each do |i|
                              puts "\t\t  #{i}     #{nbrhood_seq[i,6]}"
                         end
                    end
               end
          end
      end
end
