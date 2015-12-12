require_relative 'region'

module Utils
  ##################################
  # & is not applicable to this case, bcs two regions are equivalent if they
  # overlap, not if they have the same from and to (so which one of them should be
  # used as the representative in the intersection?)
  def find_overlapping_regions regions_list_1, regions_list_2
     overlaps = []
     uniq_for_list_1 = []
     regions_list_1.each do |region_1|
          overlap = regions_list_2.select { |region_2|  region_2  == region_1 }
          overlap.empty? ? uniq_for_list_1.push(region_1) :  overlaps.push( [region_1, overlap])
     end
     return overlaps, uniq_for_list_1
  end


  def place_pr_regions  regions, strand, exons
       
       descr_str = ""
       orig_regions = regions.dup
       if strand=='-'        
            origin  = exons[-1][1]
            regions = regions.map { |x| [-(x[1]-origin), -(x[0]-origin)]}
            regions.reverse()
            exons =  exons.map { |x| [-(x[1]-origin), -(x[0]-origin)] }
            exons.reverse()
       end
       
       first = true
       regions.each do |region|

            if strand=='-'
                 orig_index  = regions.length - regions.index(region) - 1
                 orig_region = orig_regions[orig_index]
                 descr_str += "PR binding region #{orig_region[0]} ... #{orig_region[1]}:   "
            else
                 descr_str += "PR binding region #{region[0]} ... #{region[1]}:   "
            end
            
             
            if  region[1] < exons[0][0]
                 diff = exons[0][0]-region[1]
                 if diff < 1000
                      descr_str += "#{diff} bp upstream from the first exon" 
                 else
                      descr_str += "#{format("%.1f",diff/1000.0)}  kbp upstream from the first exon"
                 end
                 
            elsif  exons[-1][1] < region[0]
                 diff = region[0] -  exons[-1][1]
                 if diff < 1000
                      descr_str += "#{diff} bp downstream from the last exon" 
                 else
                      descr_str += "#{format("%.1f",diff/1000.0)} kbp downstream from the last exon"
                 end
                 
            else
                 region_start = []
                 if  region[0] >= exons[0][0]
                      (1 ... exons.length).each do |i|
                           if region[0] <= exons[i][1]
                                region_start = ["exon", i+1]
                                break
                           end
                           if i+1 < len(exons) and region[0] < exons[i+1][0]
                                region_start = ["intron", i+1]
                                break
                           end
                      end
                 end
                 region_end = []
                 if region[0] <= exons[-1][1]
                      (1 ... exons.length).each do |i|
                           if region[1] <= exons[i][1]
                                region_end = ["exon", i+1]
                                break
                           end
                           if i+1 < exons.length and region[1] < exons[i+1][0]
                                region_end = ["intron", i+1]
                                break
                           end
                      end     
                 end
            
                 if region_start.length ==0  and region_end.length== 0
                      descr_str += "err locating PR region (?)"
                 elsif region_end.length == 0
                      descr_str +=  "starts in  #{region_start[0]}  #{region_start[1]}, ends after the last exon "
                 elsif region_start.length == 0
                      descr_str +=  "starts before the first exon, ends in #{region_end[0]}  #{region_end[1]}"
                 elsif region_start[0] == region_end[0] and  region_start[1] == region_end[1]
                      descr_str += "within  #{region_start[0]}  #{region_start[1]}" 
                 else
                      descr_str += "starts in  #{region_start[0]}  #{region_start[1]}, " 
                      descr_str += "ends in  #{region_end[0]}  #{region_end[1]}"
                 end
            end
            descr_str += "\n"
       end
       return descr_str
  end # end of place_pr_regions    
    
end # end of module
