#!/usr/bin/env ruby  -W0

require_relative 'ruby_modules/parsers'
require_relative 'ruby_modules/mysqlutils'

include Parsers, MysqlUtils


chromosomes = (1..19).map { |i| i.to_s} + ['X']

motif_positions = {}
chromosomes.each {|chrom| motif_positions[chrom] = []}
score = {}
File.readlines('best_alignments').each do |line|     
     name, scr, no_species = line.split "\s"
     chrom, start, stop = name.gsub("clean_alignments/", "").gsub(".afa","").split "_"
     motif_positions[chrom].push ([start.to_i,stop.to_i])
     id = "#{chrom}_#{start}_#{stop}"
     score[id] = scr
end

id_translation = get_id_translation # in parsers

=begin
name = "A"
puts "#{name}"
ids = id_translation[name]
ids[:other_names].each {|o| puts "\t #{o}"}
ids[:uniprots].each {|o| puts "\t #{o}"}
ids[:human_names].each {|o| puts "\t #{o}"}
ids[:human_uniprots].each {|o| puts "\t #{o}"}
exit
id_translation.to_a.sample(10).to_h.each do |name, ids|
     puts "#{name}"
     ids[:other_names].each {|o| puts "\t #{o}"}
     ids[:uniprots].each {|o| puts "\t #{o}"}
     ids[:human_names].each {|o| puts "\t #{o}"}
     ids[:human_uniprots].each {|o| puts "\t #{o}"}
end
exit
=end

# connect to ucsc, mm9 database
connection_handle = connect_to_mysql('/Users/ivana/.ucsc_mysql_conf')
connection_handle.select_db ('mm9')


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
          dist.sort_by {|k,v| v.abs}.to_h.keys[0,1].each do |name|
               next if dist[name].abs>50000
               id = "#{chrom}_#{motif_from}_#{motif_to}"
               printf "%-2s    %12d  %12d   %8.1f  %14s  %1s  %6d   ", chrom, motif_from, motif_to, score[id], name, strand[name], dist[name]
               if id_translation.has_key? name.upcase
                    ids = id_translation[name.upcase]
                    human_name = "unresolved"
                    if  ids[:human_names].length > 0 then human_name=ids[:human_names][0] end
                   
                    printf " %s  %s %s ", human_name,  ids[:uniprots].join(','), ids[:human_uniprots].join(',')
               else
                    print " unresolved "
               end
               puts  #newline
               
               if false and dist[name] == 0
                    gene_detail = {}
                    qry = "select name, name2, chrom, strand, txStart, txEnd, exonStarts, exonEnds from refGene where name2='#{name}'"
                    rows = connection_handle.query(qry)
                    if not rows or rows.size == 0 then abort "no info for #{name}" end
                    rows.each do |row|
                         gene_detail[row['name']] = {
                              chrom: row['chrom'].sub('chr',''),
                              strand: row['strand'],
                              region: Region.new(row['txStart'].to_i, row['txEnd'].to_i),
                              exonStarts: row['exonStarts'].split(',').map {|s| s.to_i},
                              exonEnds: row['exonEnds'].split(',').map {|e| e.to_i}
                         }                
                    end
                         
                    gene_detail.each do |name, data|
                         print  "\t  #{name}"
                         if data[:exonStarts][-1] < motif_from
                              puts " after the last exon"
                              break
                         end
                        (0...data[:exonStarts].length).each do |i|
                              if motif_to < data[:exonStarts][i]
                                   puts " before exon #{i} "
                                   break
                              elsif motif_to < data[:exonEnds][i]
                                   if motif_from < data[:exonStarts][i]
                                        puts " straddles the beginning of exon #{i}"
                                        break
                                   else
                                        puts " within exon #{i}"
                                        break
                                   end
                              elsif motif_from < data[:exonEnds][i]
                                   puts " straddles the end of exon #{i}"
                                   break
                              end
                        end
                    end
               end
          end
     end
end
