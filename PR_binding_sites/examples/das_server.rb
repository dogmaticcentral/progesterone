#!/usr/bin/env ruby
require 'openssl'
OpenSSL::SSL::VERIFY_PEER = OpenSSL::SSL::VERIFY_NONE
require 'rubygems'
require 'nokogiri'
require 'open-uri'
require 'open_uri_redirections'

assembly = 'mm9'
chrom = 8
from = 59791026
to = 59791040

das_request  = "http://genome.ucsc.edu/cgi-bin/das/#{assembly}/"
das_request += "dna?segment=chr#{chrom}:#{from},#{to}"

page = Nokogiri::HTML(open(das_request, :allow_redirections => :all))
puts page.class   # => Nokogiri::HTML::Document
#puts page.css("dna")[0].text
#puts
dna_seq = page.css("dna")[0].text.gsub(/[\s\n]/, '')
puts dna_seq

