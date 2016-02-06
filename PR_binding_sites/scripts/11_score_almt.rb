#!/usr/bin/env ruby -W0

ARGV.length > 0 or abort "Usage: #{$0} <afa>"

###########################################################################################
def flip_category nt
     nt=='A' and return 'G'
     nt=='G' and return 'A'
     nt=='T' and return 'C'
     nt=='C' and return 'T'
     return nt
end
# RgNaca NRN tgtNcY where R = A|G, Y = T|C, N = any
def seq_comparison_score ref_seq, other_seq
     ref_seq.length.odd? or abort "ref_seq lngth not odd"
     ref_seq.length == other_seq.length  or abort "ref_seq lngth not equal to the other seq length"
     ref_seq.length >= 19  or abort "ref_seq shorter than expected"
     score = 0
     middle_index =  ref_seq.length/2
     
     (middle_index-1 .. middle_index-1).each do |i|
          if other_seq[i] == '-'
               score -= 1
          elsif ref_seq[i] == other_seq[i]
               score += 2
          elsif ref_seq[i]==flip_category(other_seq[i])
               score += 1
          end
     end
     (2 .. 4).each do |offset|
          
          [middle_index + offset, middle_index-offset].each do |i|
           
               if other_seq[i] == '-'
                    score -= 4
               elsif ref_seq[i] == other_seq[i]
                    score += 4
               elsif ref_seq[i]==flip_category(other_seq[i])
                    score += 1
               end
          end
     end
     offset = 5
     [middle_index + offset, middle_index-offset].each do |i|
          if other_seq[i] == '-'
               score -= 1
          elsif ref_seq[i] == other_seq[i]
               score += 1
          elsif ref_seq[i]==flip_category(other_seq[i])
               score += 0.5
          end
     end
     offset = 6
     [middle_index + offset, middle_index-offset].each do |i|
          if other_seq[i] == '-'
               score -= 3
          elsif ref_seq[i] == other_seq[i]
               score += 3
          elsif ref_seq[i]==flip_category(other_seq[i])
               score += 2
          end
     end
     offset = 7
     [middle_index + offset, middle_index-offset].each do |i|
          if other_seq[i] == '-'
               score -= 2
          elsif ref_seq[i] == other_seq[i]
               score += 2
          elsif ref_seq[i]==flip_category(other_seq[i])
               score += 1
          end
     end
     (8 .. middle_index).each do |offset|
          [middle_index + offset, middle_index-offset].each do |i|              
              score += 0.5 if ref_seq[i] == other_seq[i] 
          end
     end
          
     return score
         
end
###########################################################################################
def remove_gaps! alignment
     # remove gapped positions
     k,v = alignment.first
     almt_length = v.length
     if not almt_length > 0 then return end # not my job
     gapped = []
     (0 ... almt_length).each  do |i|
          gapped_in = alignment.select {|k,v| v[i] == '-'}
          next if gapped_in.length != alignment.length
          gapped.push(i)
     end
     if gapped.length == 0 then return end
     alignment.each_key do |k|
          new_seq = ""
          (0 ... almt_length).each { |i| new_seq += alignment[k][i] if not gapped.include? i}
          alignment[k] = new_seq          
     end
end
###########################################################################################
def cleanup! alignment, ref_species, relevant_species
     # stick to interesting species only
     alignment.keep_if {|k,v| k==ref_species or relevant_species.include? k}
     maxlength = alignment.values.map {|v| v.length}.max
     # not sure if I understant this corretly, but it looks like whichever thing is giving me alignments
     # does not pad to the full sequence length if there is a series of gaps
     alignment.each_key  {|k|  alignment[k] += '-'*(maxlength - alignment[k].length)}
     remove_gaps! alignment
     # if the ref sequence contains no gaps, we are done
     alignment[ref_species].include? '-' or return
     # remove species with inserts wrt to ref qry
     ref_seq = alignment[ref_species]
     alignment.keep_if do |k,v|
          ( (0 ... ref_seq.length).each.select {|i| ref_seq[i]=='-' and v[i] != '-'}.length) == 0
     end
     remove_gaps! alignment
     return    
end
###########################################################################################
require_relative 'ruby_modules/httputils'
include  HttpUtils


file_name = ARGV[0]

seq = {}
name = ""
File.readlines(file_name).each do |line|
     line.chomp!
     if line[0]=='>'
          name, chrom, from, to = line[1..-1].split ' '
          seq[name] = ""
     else
          seq[name] += line
     end
end

species = [ "anoCar1", "bosTau3", "calJac1", "canFam2", "cavPor2",
           "danRer5", "dasNov1", "echTel1", "equCab1", "eriEur1",
           "felCat3", "fr2", "galGal3", "gasAcu1", "hg18", "loxAfr1",
           "mm9", "monDom4", "ornAna1", "oryCun1", "oryLat1", "otoGar1",
           "panTro2", "ponAbe2", "rheMac2", "rn4", "sorAra1",
           "tetNig1", "tupBel1", "xenTro2"]
# is this the right seq of species/assembly names?
seq.each  do |seqname, sequence|
     species.include? seqname or abort "Unrecognized species (or assembly) #{seqname}"
end


#  how do we want to group the species
ref_species = "mm9"
species_names = {
     rodents:  ["rn4", "cavPor2", "oryCun1"],
     primates: ["hg18", "panTro2", "ponAbe2", "rheMac2",  "calJac1"],
     other_placentals: [ "bosTau3", "canFam2","dasNov1", "echTel1", "equCab1",
                  "eriEur1","felCat3","loxAfr1",  "otoGar1","sorAra1","tupBel1"],
     other_mammals: ["ornAna1","monDom4"],
     other_verts: [ "anoCar1", "danRer5", "galGal3", "gasAcu1", "fr2","oryLat1","tetNig1", "xenTro2"]
}

# cleanup the alignment: remove species that we will not be considering and realign
# what if I still have gaps? remove species that do not have the gaps in the same place
# as the reference sequence, and realign
cleanup! seq, ref_species, species_names[:rodents]+species_names[:primates]+species_names[:other_placentals]
if seq.length>1
     fnm = file_name.split('/').pop
     outf = open("clean_alignments/#{fnm}", 'w')
     seq.each {|k,v| outf.write ">#{k}\n#{v}\n"}
     outf.close
end

# keep only the species that  we have present in the alignment
[:rodents, :primates, :other_placentals].each do |species_group|
     species_names[species_group].select! {|s| seq.keys.include? s}
end

# score the alignment according to the weight we give to individual positions in
score = {}
[:rodents, :primates, :other_placentals].each do |species_group|
     score[species_group] = 0
     species_names[species_group].each do |spec|
          score[species_group] += seq_comparison_score seq[ref_species], seq[spec]
          #printf "\t\t  %-10s   %6.1f \n", spec,  seq_comparison_score(seq[ref_species], seq[spec])
     end
     species_names[species_group].length > 0 or next
     score[species_group] /= 1.0*species_names[species_group].length
     #printf "  %-10s   %6.1f \n", species_group, score[species_group]
end


total_score = 3*score[:rodents] + 2*score[:primates] + score[:other_placentals]

printf " %30s    %6.1f   %3d  \n", file_name, total_score, seq.length
=begin
 Erinaceus europaeus = European hedgehog
 Sorex araneus = Common shrew
 glires: rodents and rabbits
 Oryctolagus cuniculus  = European Rabbit
 Echinops telfairi (small Madagascar hedgehog)
Ornithorhynchus anatinus = platypus
Dasypus novemcinctus (nine-banded armadillo)
Tupaia belangeri (northern tree shrew)

=end
