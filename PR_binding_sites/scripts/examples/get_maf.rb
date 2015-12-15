#!/usr/bin/env ruby
require 'openssl'
OpenSSL::SSL::VERIFY_PEER = OpenSSL::SSL::VERIFY_NONE
require 'rubygems'
require 'nokogiri'
require 'open-uri'
require 'open_uri_redirections'

assembly = 'mm9'
chrom = 8
from = 59798898
to = 59798903

request  = "https://genome.ucsc.edu/cgi-bin/hgTables?clade=mammal&org=Mouse&"
request += "db=#{assembly}&hgta_group=compGeno&hgta_track=cons30way&hgta_table=multiz30way&gta_regionType=userRegions&"
request += "hposition=chr#{chrom}%:#{from}-#{to}&"
request += "hgta_outputType=wigData"
maf = Nokogiri::HTML(open(request, :allow_redirections => :all))

puts maf
