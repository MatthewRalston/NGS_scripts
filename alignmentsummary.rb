#!/home/mrals/.rvm/rubies/ruby-2.0.0-p247/bin/ruby
=begin
------------------------------------------------------------------------------
--                                                                          --
--                                 MATTEW RALSTON                           --
--                                                                          --
--                               S K E L E T O N . R B                      --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Summer 2013                                                             --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed ...
--                                                                          --
------------------------------------------------------------------------------
=end

################################################
#
#               R E Q U I R E
#
################################################





################################################
#
#               U S E R    V A R I A B L E S
#
################################################





################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################



def main
  x=0; rrnarate=[]; postfilterrate=[]; total=[]; minusrrna=[]
  begin
    puts('hello')
    while line=STDIN.gets
      if x%2 == 0 and (line.include?('reads; of these:') or line.include?('overall alignment rate'))
        total << line.split[0] if line.include?('reads; of these:')
        (rrnarate << line.split('%')[0]; x+=1) if line.include?('overall alignment rate')
      elsif x%2 == 1 and (line.include?('reads; of these:') or line.include?('overall alignment rate')) 
        minusrrna << line.split[0] if line.include?('reads; of these:')
        (postfilterrate << line.split('%')[0]; x+=1) if line.include?('overall alignment rate')
      end
    end
  rescue NoMethodError
  end
  

  rrnarate.size.times do |c|
    STDOUT.puts([rrnarate[c], postfilterrate[c], total[c], minusrrna[c]].join("\t"))
  end
end
  

#*****************************************************************************#
################################################
#
#-----------------------------------------------
#
#                   M A I N
#-----------------------------------------------
#
################################################



main



##########################  E O F   ############################################
