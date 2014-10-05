#!/home/mrals/.rvm/rubies/ruby-2.0.0-p247/bin/ruby
=begin
------------------------------------------------------------------------------
--                                                                          --
--                                MATTHEW RALSTON                           --
--                                                                          --
--                              C O N T A I N E D. R B                      --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Summer 2013                                                             --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed to identify items in a list that are contained in  --
-- other items of the list. This is meant to find the transcripts that are  --
-- fully contained in other transcripts.
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

GTF="Trinity_ref/Trinity_filtered.gtf"
OUTPUT="Trinity_ref/Trinitycurated.gtf"


################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################


def gtfread(infile)
  liszt=[]
  while line=STDIN.gets 
    liszt << line.chomp.split
  end
  liszt.map! {|x| x[3]=x[3].to_i; x[4]=x[4].to_i; x}
  liszt.sort! {|x,y| x[3]<=>y[3]}
end

def contained?(l1,l2)
  contained=[]
  l1.each do |g1|
    l2.each do |g2|
      if g2[3] <= g1[3] && g1[4] <= g2[4] && g1[0] == g2[0] && g1[-1] == g2[-1] && g1 != g2
        contained << g1
        break
      end
    end
  end
  contained
end

def gtfout(records,liszt,output)
  records.delete_if {|x| liszt.include?(x)}
  records.each do |x|
    STDOUT.puts(x.join("\t"))
  end
end

def main
  gtf=gtfread(GTF)
  contained=contained?(gtf,gtf)
  gtfout(gtf,contained,OUTPUT)
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
