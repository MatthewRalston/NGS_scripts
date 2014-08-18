#!/home/mrals/.rvm/rubies/ruby-2.0.0-p247/bin/ruby
=begin
------------------------------------------------------------------------------
--                                                                          --
--                                MATTHEW RALSTON                           --
--                                                                          --
--                         C O U N T S U M M A R Y . R B                    --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Copyright 2014 Matthew Ralston                                          --
--                                                                          --
------------------------------------------------------------------------------
-- This script is designed to summarize the raw read counts per gene        --
-- produced from HTSeq-count. This file condenses them into a tab delimited --
-- count matrix printed to countsummary.txt				    --
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

DIR="/home/mrals/Final/Expression/counts"



################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################


def parse(file,dict)
  File.open(file,'r').each do |line|
    if line.include?("CA_")
      dict[line.split[0]] ? dict[line.split[0]]+=line.chomp.split[1].to_i : dict[line.split[0]]=line.chomp.split[1].to_i
    end
  end
  return dict
end


def main
  dict={};new={}
  Dir.foreach(DIR) do |file|
    dict[file.split(".")[0]] = {} unless dict[file.split(".")[0]]
    dict[file.split(".")[0]] = parse(DIR+"/"+file,dict[file.split(".")[0]]) if file.include?(".counts")
  end
  conds=[]
  dict.each do |cond,dict2|
    conds << cond
    dict2.each do |gene,count|
      new[gene]={} unless new[gene]
      new[gene][cond]=count 
    end
  end
  conds.shift
  new.each do |gene, counts|
    conds.each do |cond|
      new[gene][cond]="N/A" unless new[gene][cond]
    end
  end

  File.open("summary/countsummary.txt",'w') do |file|
    file.puts("Gene_id\ttreatment\ttime\treplicate\tcounts")
    new.each do |gene,counts|
      counts.each do |cond,count|
        file.puts("#{gene}\t#{cond}\t#{cond}\t#{cond}\t#{count}")
      end
    end
  end
  new["CA_C0001"].keys.each do |condition|
    File.open("Expression/"+condition+".3.counts",'w') do |file|
      new.each do |gene,counts|
        file.puts("#{gene}\t#{counts[condition]}")
      end
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
