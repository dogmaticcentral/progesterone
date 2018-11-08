#!/usr/bin/env ruby
require 'openssl'
OpenSSL::SSL::VERIFY_PEER = OpenSSL::SSL::VERIFY_NONE
require 'rubygems'
require 'nokogiri'
require 'open-uri'
require 'open_uri_redirections'

module HttpUtils
     def get_dna_region assembly, chrom, from, to
          das_request  = "http://genome.ucsc.edu/cgi-bin/das/#{assembly}/"
          das_request += "dna?segment=chr#{chrom}:#{from},#{to}"
          page = Nokogiri::HTML(open(das_request, :allow_redirections => :all))
          dna_seq =  page.css("dna")[0].text.gsub(/[\s\n]/, '')
          return dna_seq
     end

     def get_maf assembly, chrom, from, to
          request  = "https://genome.ucsc.edu/cgi-bin/hgTables?clade=mammal&org=Mouse&"
          request += "db=#{assembly}&hgta_group=compGeno&hgta_track=cons30way&hgta_table=multiz30way&gta_regionType=userRegions&"
          request += "hposition=chr#{chrom}%:#{from}-#{to}&"
          request += "hgta_outputType=wigData"
          maf = Nokogiri::HTML(open(request, :allow_redirections => :all))
          return maf
     end
end

