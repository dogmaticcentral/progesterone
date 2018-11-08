#!/usr/bin/env ruby

include Math

################################
def parse_ret ret, qry_from, qry_to
     motif = nil
     motif_score = nil
     distance = nil
     strand = nil
     chrom = nil
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
               motif_score = score
          end
     end
     return motif, motif_score.to_f, distance.to_i, strand, chrom
end
################################
def expression idx, gene_name
     if idx ==1
          table = "/Users/ivana/Dropbox/Biserka/PR_ChIPSeq_paper/SM3_table_1:Microarray_run1.csv"
     else
          table = "/Users/ivana/Dropbox/Biserka/PR_ChIPSeq_paper/SM3_table_2:Microarray_run2.csv"
     end
     ret = `grep #{gene_name} #{table}`
     if ret.length == 0
          return '-'
     else
          if idx==1
               return sprintf "%5.2f", ret.split(' ')[2].to_f
          else
               return sprintf "%5.2f", ret.split(' ')[1].to_f
          end
     end
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

score  = score0.sort_by {|k,v|  -v}

max_alnmt_score = -1
max_motif_score = -1
min_dist        =  10000000
max_dist        = -1
data = {}
score.each do |gene_name, alnmt_score|
     motif_op = nil
     motif_p4 = nil
     motif    = nil
     cond  = nil
     chrom = 0
     motif_score = 0
     ret_op =  `grep #{gene_name} top_genes_half_site.oil_and_p4.stricter`
     ret_p4 =  `grep #{gene_name} top_genes_half_site.p4_only.stricter`
    
     if ret_op.length>0
          motif_op, motif_score, distance, strand, chrom = parse_ret(ret_op, from[gene_name], to[gene_name])
          cond  = "o&p"         
     end
     if  ret_p4.length > 0
          motif_p4, motif_score, distance, strand, chrom = parse_ret(ret_p4, from[gene_name], to[gene_name])
          cond = "p4 "
     end

     if motif_op and motif_p4
          abort "\nsame motif in o&p and p4 only(?)"
     elsif motif_op
          motif = motif_op         
     elsif motif_p4
          motif = motif_p4
     else
          abort "\nno info found in top_genes table (?)"
     end
     
     data[gene_name] = [alnmt_score, seqs_in_almt[gene_name], motif, motif_score, cond, distance, strand, chrom]
     max_alnmt_score < alnmt_score and max_alnmt_score = alnmt_score
     max_motif_score < motif_score and max_motif_score = motif_score
     min_dist > distance and min_dist = distance
     max_dist < distance and max_dist = distance
end

cumulative_score = {}
data.each do |gene_name, d|
     alnmt_score, no_of_seqs, motif, motif_score, cond, distance = d
     cs  = ( (max_alnmt_score-alnmt_score)/max_alnmt_score )**2
     cs += ( (max_motif_score-motif_score)/max_motif_score )**2
     cs += ( (distance-min_dist)/max_dist )**2
     cumulative_score[gene_name]  = sqrt(cs)
end

gene_names_sorted =  data.keys.sort_by {|v|  cumulative_score[v]}

gene_names_sorted.each do |gene_name|
     
     alnmt_score, no_of_seqs, motif, motif_score, cond, distance, strand, chrom = data[gene_name]
     printf " %15s   %5.2f  %3d ", gene_name, alnmt_score, seqs_in_almt[gene_name]
     printf " %3s %s %18s  %5.2f  %3s   %6d   %5.3f  ", chrom, strand, motif, motif_score,
     cond, distance, cumulative_score[gene_name]
     printf "   %5s     %5s  \n", expression(1, gene_name),  expression(2, gene_name)
end
