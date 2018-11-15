#!/usr/bin/env ruby

require_relative 'ruby_modules/mysqlutils'

include MysqlUtils

connection_handle = connect_to_mysql('/Users/ivana/.ucsc_myql_conf')

connection_handle.select_db ('mm9')

qry = 'show tables like "ref%"'
puts "="*20
puts qry
ret = connection_handle.query qry
ret.each do |row|
     row.each { |k,v| puts " #{k}: #{v}" }
end

qry = 'select * from refGene limit 10'
puts "="*20
puts qry
connection_handle.query(qry).each do |row|
     row.each { |k,v| puts " #{k}: #{v}" }
end

connection_handle.close
