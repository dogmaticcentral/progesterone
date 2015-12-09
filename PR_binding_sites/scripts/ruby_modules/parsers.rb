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
     def parse_gene_table table_name, extension=10000
          gene_regions = []
          gene_data    = {}
          File.readlines(table_name).drop(1).each do |line|
               name, name2, txStart, txEnd = line.split "\t"
               region = Region.new(txStart.to_i-extension,txEnd.to_i+extension)
               gene_regions.push(region)
               gene_data[region] = line
          end
          return gene_data, gene_regions
     end
end
