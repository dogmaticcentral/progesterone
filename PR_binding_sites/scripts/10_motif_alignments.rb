#!/usr/bin/env ruby

require_relative 'ruby_modules/region'
require_relative 'ruby_modules/httputils'
include  HttpUtils

#table_name                 =  "top_genes_half_site.p4_only.stricter"
table_name                 =  "top_hand2"
maf_region_extraction_tool =  "/Users/ivana/third_party_utils/kentUtils/bin/mafsInRegion"
maf_dir                    =  "/Users/ivana/databases/UCSC"
maf2afa_tool               =  "/usr/local/bin/maf_to_fasta.py"
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

     # now get this region from the alignment
     outf = open("regions.bed", 'w')
     outf.write("chr#{gene.chrom}  #{motif_from-1}  #{motif_to}\n") 
     outf.close
     system ("#{maf_region_extraction_tool}  regions.bed   tmp.maf   #{maf_dir}/chr#{gene.chrom}.maf")
     
     if not File.exist? "tmp.maf"
          puts "maf not produced"
          next
     end
     #turn to afa
     system ("#{maf2afa_tool} < tmp.maf > tmp.afa")
     if not File.exist? "tmp.afa"
          puts "afa not produced"
          next
     end
     if  File.zero? "tmp.afa"
          puts "afa empty"
          next
     end     
     system ("rm tmp.maf")

     
     
     # read in, stich the pieces of the sequnce, because of the idiotic format in which the script returns it
     aligned_seq = {}
     last_pos = {}
     first_pos = {}
     key = ""
     File.readlines("tmp.afa").each do |line|
          if line[0]=='>'
               crap, from_to = line[1..-1].split (/\:/)
               from, to = from_to.split /\-/
               fields = crap.split /\./
               assembly = fields.shift
               chrom = fields.join '.'
               
               key = "#{assembly} #{chrom}"
               if not aligned_seq.keys.include? key 
                    aligned_seq[key] = ""
                    first_pos[key] = to.to_i                    
               elsif from.to_i != last_pos[key]
                     aligned_seq[key] += "_pos_mismatch_: #{from} #{last_pos[key]}   #{key}"
               end
               last_pos[key] = to.to_i
          else
               if key=="" then abort "key not initialized" end
               aligned_seq[key] += line.chomp
          end
     end
     system ("rm tmp.afa")
    puts

     
     outf = open(alignment_file, 'w')
     aligned_seq.each do |seqname, sequence|
          # cleanup
          next if sequence.include? "_pos_mismatch_"
          next if sequence.gsub('-','').length < query_seq.length
          outf.write(">#{seqname} #{first_pos[seqname]} #{last_pos[seqname]}\n")
          outf.write(sequence+"\n")
     end
     outf.close
     # extract rodents and primates  and ungulates (?) if available
     
end

