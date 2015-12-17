#!/usr/bin/env ruby

def parse_ret ret, qry_from, qry_to
     motif = nil
     alignment_score = nil
     distance = nil
     ret.split("\n").each do |line|
          name, strand, chrom, f, t, distance, length, index, motif_l, middle, right, score = line.split
          if  strand == '+'
               motif_from = f.to_i  - distance.to_i  - length.to_i  + index.to_i - 9
          else
               motif_from = t.to_i  + distance.to_i +  index.to_i - 9
          end
          motif_to   = motif_from + 18
          if motif_from == qry_from  and motif_to == qry_to
               motif and abort "motif found already (?)"
               motif = "#{motif_l} #{middle} #{right}"
               alignment_score = score
          end
     end
     return motif, alignment_score, distance
end

################################

score0 = {}
seqs_in_almt = {}
from = {}
to = {}
File.readlines('almt_score').each do |line|
     line.chomp!
     name, sc, number_of_seqs = line.split
     next if  sc.to_f <220 or number_of_seqs.to_i <12
     name_fields =  name.sub('alignments/','').sub('.afa','').split('_')
     gene_name  = name_fields[0]
     from[gene_name] = name_fields[1].to_i
     to[gene_name]   = name_fields[2].to_i
     score0[gene_name] = sc.to_f
     seqs_in_almt[gene_name] = number_of_seqs.to_i
end

score = score0.sort_by {|k,v|  -v}

score.each do |gene_name,v|
     printf " %15s   %5.2f  %3d ", gene_name, v,  seqs_in_almt[gene_name]
     motif_op = nil
     motif_p4 = nil
     cond  = nil
     ret_op =  `grep #{gene_name} top_genes_half_site.oil_and_p4.stricter`
     ret_p4 =  `grep #{gene_name} top_genes_half_site.p4_only.stricter`

     
     if ret_op.length>0
          motif_op, almt_score_op, distance = parse_ret(ret_op, from[gene_name], to[gene_name])
          cond  = "o&p"
     end
     if  ret_p4.length > 0
          motif_p4, almt_score_p4, distance = parse_ret(ret_p4, from[gene_name], to[gene_name])
          cond = "p4 "
     end

     if motif_op and motif_p4
          abort "\nsame motif in o&p and p4 only(?)"
     elsif motif_op
          printf " %20s  %5.2f  %3s   %6d \n", motif_op, almt_score_op, cond, distance
          
     elsif motif_p4
          printf " %20s  %5.2f  %3s   %6d \n", motif_p4, almt_score_p4, cond, distance
     else
          abort "\nno info found in top_genes table (?)"
     end

     
end
