##################################
class Region
     attr_accessor :from, :to
     def initialize from, to
          @from = from
          @to = to
     end

     def  overlaps? reg2
          not (@to < reg2.from or reg2.to < @from)
     end

     def covers? reg2
          @from < reg2.from and reg2.to < @to
     end

     def is_sub? reg2
          reg2.from < @from and @to < reg2.to
     end

     def is_identical? reg2
          reg2.from == @from and @to == reg2.to
     end

     alias === is_identical?
     alias == overlaps?
     alias > covers?
     alias < is_sub?

     def merge reg2
          new_from = @from<reg2.from ? @from : reg2
          new_to   = @to>reg2.to ? @to: reg2.to
          return Region.new(new_from, new_to)
     end

     def length
          return @to - @from + 1
     end
end

