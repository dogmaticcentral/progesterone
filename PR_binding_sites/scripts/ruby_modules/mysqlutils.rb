require 'mysql'

module MysqlUtils

     def connect_to_mysql (conf_file)
          # a primitive parser that looks for two or three fields that it expects
          host = nil 
          user = nil
          flag = nil
          File.readlines(conf_file).each  do |line|
               if line.include? "host" then  host  = line.split.pop.gsub '"',''
               elsif  line.include? "user" then user = line.split.pop.gsub '"',''
               elsif line.include? "skip-auto-rehash" then flag = "-A"
               end
          end
          puts host, user, flag

          conn_handle = nil
          begin
               conn_handle = Mysql.new host=host, user=user, flag=flag
               puts conn_handle.get_server_info
               rs = conn_handle.query 'SELECT VERSION()'
               puts rs.fetch_row    
               
          rescue Mysql::Error => e
               puts e.errno
               puts e.error
    
          ensure
               conn_handle.close if conn_handle
               conn_handle = nil 
          end

     end
     
end
