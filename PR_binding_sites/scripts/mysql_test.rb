#!/usr/bin/env ruby

require_relative 'ruby_modules/mysqlutils'

include MysqlUtils

connection_handle = connect_to_mysql('/Users/ivana/.ucsc_myql_conf')
puts connection_handle
