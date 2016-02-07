#!/usr/bin/env ruby  -W0

require_relative 'ruby_modules/parsers'
require_relative 'ruby_modules/mysqlutils'

include Parsers, MysqlUtils


chromosomes = (1..19).map { |i| i.to_s} + ['X']

motif_positions = {}
chromosomes.each {|chrom| motif_positions[chrom] = []}
     
File.readlines('best_alignments').each do |line|
     
     name, score, no_species = line.split "\s"
     chrom, start, stop = name.gsub("alignments/", "").gsub(".afa","").split "_"
     motif_positions[chrom].push ([start.to_i,stop.to_i])
end
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
          dist.sort_by {|k,v| v.abs}.to_h.keys[0,3].each do |name|
               next if dist[name].abs>50000
               puts "#{chrom}   #{motif_from} #{motif_to}   #{name}  #{strand[name]}    #{dist[name]}"
               
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
