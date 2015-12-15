#!/usr/bin/env ruby
require 'openssl'
OpenSSL::SSL::VERIFY_PEER = OpenSSL::SSL::VERIFY_NONE
require 'rubygems'
require 'nokogiri'
require 'open-uri'
require 'open_uri_redirections'


das_request  = "http://genome.ucsc.edu/cgi-bin/das/hg19/"
das_request += "dna?segment=chr1:100000,100500"

page = Nokogiri::HTML(open(das_request, :allow_redirections => :all))
puts page.class   # => Nokogiri::HTML::Document
puts page.css("dna")[0].text
puts
dna_seq = page.css("dna")[0].text.gsub(/[\s\n]/, '')
puts dna_seq

