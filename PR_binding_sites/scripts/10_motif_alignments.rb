#!/usr/bin/env ruby

require_relative 'ruby_modules/region'
require_relative 'ruby_modules/httputils'
include  HttpUtils

#table_name                 =  "top_genes_half_site.p4_only.stricter"
table_name                 =  "top_hand2"
class Gene
     attr_accessor :stradn, :chrom, :region, :strand, :chrom, :motif_index
     attr_accessor :dist_to_PR_binding_region, :length_of_PR_binding_region, :motif, :score
     def initialize strand, chrom, from, to, dst_to_PR_binding, binding_length, motif_index, motif_left, middle, right, score
          @strand = strand
          @chrom  = chrom
          @region = Region.new(from.to_i,to.to_i)
          @dist_to_PR_binding_region   = dst_to_PR_binding.to_i
          @length_of_PR_binding_region = binding_length.to_i
          @motif_index = motif_index.to_i
          @motif =  "#{motif_left} #{middle} #{right}"
          @score = score
     end

end
#chr11:68910576-68910594

genes = {}
File.readlines(table_name).each do |line|
     gene_name, strand, chrom, from, to, dst_to_PR_binding, binding_length, motif_index, motif_left, middle, right, score = line.split
     genes[gene_name] = Gene.new strand, chrom, from, to, dst_to_PR_binding, binding_length, motif_index, motif_left, middle, right, score
     
end
reversed_genes = Hash[genes.to_a.reverse]

genes.each do |gene_name, gene|
#if true
     #gene_name = 'Scnn1a'
     #gene = genes[gene_name]
     puts " #{gene_name}   #{gene.chrom}  #{gene.score}"
     #checking
     if  gene.strand == '+'
          motif_from = gene.region.from - gene.dist_to_PR_binding_region - gene.length_of_PR_binding_region + gene.motif_index - 9
     else
          motif_from = gene.region.to  +  gene.dist_to_PR_binding_region + gene.motif_index - 9
     end
     motif_to   = motif_from + 18
     
     alignment_file = "alignments/#{gene_name}_#{motif_from}_#{motif_to}.afa"
     next if File.exist? alignment_file
  
     seq =  get_dna_region 'mm9', gene.chrom, motif_from, motif_to
     query_seq = gene.motif.gsub(" ", "")
     puts seq
     puts "  "+query_seq
     if seq[2..-3] != query_seq then  abort "seq mismatch" end
     # get_alignment is in utils
     get_alignment alignment_file, gene.chrom  motif_from-1,   motif_to

     
end

