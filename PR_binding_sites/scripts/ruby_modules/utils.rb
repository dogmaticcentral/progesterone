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
end
