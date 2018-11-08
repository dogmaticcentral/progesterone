#!/usr/bin/env ruby
require 'openssl'
OpenSSL::SSL::VERIFY_PEER = OpenSSL::SSL::VERIFY_NONE
require 'rubygems'
require 'nokogiri'
require 'open-uri'
require 'open_uri_redirections'

page = Nokogiri::HTML(open("https://en.wikipedia.org/", :allow_redirections => :all))
puts page.class   # => Nokogiri::HTML::Document
puts page.css("title")[0].text

