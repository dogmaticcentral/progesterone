require 'mysql2'

module MysqlUtils

     def connect_to_mysql (conf_file)
          conn_handle = nil
          begin
               conn_handle = Mysql2::Client.new(:default_file => conf_file)
               rs = conn_handle.query 'SELECT VERSION()'
               rs.each do |row|
                    row.each { |k,v| puts "connected to #{k}: #{v}" }
               end
 
          rescue Mysql2::Error => e
               puts e.errno
               puts e.error
    
               conn_handle.close if conn_handle
               conn_handle = nil 
          end
          return conn_handle
     end
     
end
