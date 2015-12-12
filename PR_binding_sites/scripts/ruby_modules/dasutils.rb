#!/usr/bin/env ruby
require 'openssl'
OpenSSL::SSL::VERIFY_PEER = OpenSSL::SSL::VERIFY_NONE
require 'rubygems'
require 'nokogiri'
require 'open-uri'
require 'open_uri_redirections'

module DasUtils
     def get_dna_region assembly, chrom, from, to
          das_request  = "http://genome.ucsc.edu/cgi-bin/das/#{assembly}/"
          das_request += "dna?segment=chr#{chrom}:#{from},#{to}"

          page = Nokogiri::HTML(open(das_request, :allow_redirections => :all))
          dna_seq =  page.css("dna")[0].text.gsub(/[\s\n]/, '')
          return dna_seq
     end
end

