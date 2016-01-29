require_relative 'region'

module Parsers
     ##################################
     def parse_chipseq_table table_name
          chromosome = Hash.new
          File.readlines(table_name).drop(1).each do |line|
               sample, id, chrom, from, to = line.split "\t"
               (chromosome.key? chrom)  || chromosome[chrom] = Array.new 
               chromosome[chrom].push(Region.new(from.to_i,to.to_i))
          end
          return chromosome
     end
     
     ##################################
     def parse_gene_table(table_name: nil, extension: 0, merge_splices: false)

          if not table_name then abort "no table name given in parse_gene_table()" end
          gene_regions = []
          gene_data    = {}
          raw_input    = {}
          File.readlines(table_name).drop(1).each do |line|
               name, name2, strand, txStart, txEnd = line.split "\t"
               if merge_splices
                    raw_input.keys.include? name2 or raw_input[name2] = []
                    raw_input[name2].push ( [name, strand, txStart, txEnd] )
               else
                    region = Region.new(txStart.to_i-extension,txEnd.to_i+extension)
                    gene_regions.push(region)
                    gene_data[region] = line
               end
          end
          if merge_splices
               raw_input.each do |name, data|
                    splice_names = data.map {|d| d[0]}.join ","
                    strands      = data.map {|d| d[1]}.uniq
                    if strands.length>1 then
                         #puts "  #{table_name}  "
                         #puts  "strands not uniq (?)   #{name}  #{strands.inspect}"
                         next
                    end
                    strand  = strands[0]
                    tx_from = data.map {|d| d[2].to_i}.min
                    tx_to   = data.map {|d| d[3].to_i}.max
                    
                    region = Region.new(tx_from-extension,tx_to+extension)
                    gene_regions.push(region)
                    gene_data[region] = [name, splice_names, strand, tx_from, tx_to]
               end
          end
          return gene_data, gene_regions
     end
     
     ####################################
     def read_dna condition, chrom
          file = File.open("../data_raw/#{condition}_chr#{chrom}.fasta", "r")
          dna_seqs = {}
          name = ""
          file.each_line do |line|
               line.chomp!.strip!
               next if line.length == 0
               if line[0] == '>'
                    name = line[1..-1].chomp
                    dna_seqs[name] = ""
               else
                    dna_seqs[name] += line
               end
          end
          file.close
          return dna_seqs
     end

     
end
