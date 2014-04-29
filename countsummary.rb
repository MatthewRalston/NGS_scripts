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

DIR="/home/mrals/ETP/counts"



################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################


def parse(file)
  dict={}
  File.open(file,'r').each do |line|
    dict[line.split[0]]=line.split[1] if line.include?("CA_")
  end
  return dict
end


def main
  dict={};new={}
  Dir.foreach(DIR) do |file|
    dict[file[0...-9]] = parse(DIR+"/"+file) if file.include?(".counts")
  end
  conds=[]
  dict.each do |cond,dict2|
    conds << cond
    dict2.each do |gene,count|
      new[gene]={} unless new[gene]
      new[gene][cond]=count 
    end
  end
  new.each do |gene, counts|
    conds.each do |cond|
      new[gene][cond]="N/A" unless new[gene][cond]
    end
  end

  File.open("countsummary.txt",'w') do |file|
    file.puts("Gene_id\t"+new["CA_C0001"].keys.join("\t"))
    new.each do |gene,counts|
      file.puts(([gene]+counts.values).join("\t"))
    end
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
