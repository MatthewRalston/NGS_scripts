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


def genome(liszt)
=begin
NSblahblahblah
7621912 reads; of these:
  4009241 (52.60%) were paired; of these:
    1983557 (49.47%) aligned concordantly 0 times
    1818685 (45.36%) aligned concordantly exactly 1 time
    206999 (5.16%) aligned concordantly >1 times
    ----
    1983557 pairs aligned concordantly 0 times; of these:
      409503 (20.64%) aligned discordantly 1 time
    ----
    1574054 pairs aligned 0 times concordantly or discordantly; of these:
      3148108 mates make up the pairs; of these:
        1009275 (32.06%) aligned 0 times
        35392 (1.12%) aligned exactly 1 time
        2103441 (66.82%) aligned >1 times
  3612671 (47.40%) were unpaired; of these:
    498719 (13.80%) aligned 0 times
    2246121 (62.17%) aligned exactly 1 time
    867831 (24.02%) aligned >1 times
=end
  dict={}; liszt.shift
  dict["total"]=liszt.shift.split[0]; liszt.shift
  dict["paired"]=liszt.shift.split[0]
  dict["conc_once"]=liszt.shift.split[0]
  dict["conc_mult"]=liszt.shift.split[0]
  liszt.shift(2); dict["disc_once"]=""; dict["disc_mult"]=""
  line=liszt.shift
  line.include?(">1 times") ? dict["disc_mult"]=line.split[0] : dict["disc_once"]=line.split[0]
  liszt.shift
  dict["unaligned_pairs"]=liszt.shift.split[0]
  liszt.shift
  dict["unmates"]=liszt.shift.split[0]
  dict["mate_once"]=liszt.shift.split[0]
  dict["mate_mult"]=liszt.shift.split[0]
  dict["unpaired"]=liszt.shift.split[0]
  dict["unpair_unaligned"]=liszt.shift.split[0]
  dict["unpair_once"]=liszt.shift.split[0]
  dict["unpair_mult"]=liszt.shift.split[0]
  dict
end

def pair(liszt)
=begin
NSblahblahblah
12931751 reads; of these:
  12931751 (100.00%) were paired; of these:
    4009241 (31.00%) aligned concordantly 0 times
    5093 (0.04%) aligned concordantly exactly 1 time
    8917417 (68.96%) aligned concordantly >1 times
=end
  dict={}; liszt.shift
  dict["total"]=liszt.shift.split[0]; liszt.shift
  dict["unaligned"]=liszt.shift.split[0]
  dict["once"]=liszt.shift.split[0]
  dict["multiple"]=liszt.shift.split[0]
  dict
end

def unpair(liszt)
=begin
NSblahblahblah
12377949 reads; of these:
  12377949 (100.00%) were unpaired; of these:
    3612671 (29.19%) aligned 0 times
    12286 (0.10%) aligned exactly 1 time
    8752992 (70.71%) aligned >1 times
=end
  dict={}; liszt.shift
  dict["total"]=liszt.shift.split[0]
  dict["unaligned"]=liszt.shift.split[0]
  dict["once"]=liszt.shift.split[0]
  dict["multiple"]=liszt.shift.split[0]
  dict
end

def printer(dict,name)
  File.open(name,'w') do |file|
    file.puts("ID\t"+dict[dict.keys[0]].keys.join("\t"))
    dict.each {|key,value| file.puts(key+"\t"+value.values.join("\t"))}
  end
end

def main
  x=0; pairrrnarate=[];unrrnarate=[]; alignment=[]
  temp=[]; total={};paired={};unpaired={}
  begin

    while line=STDIN.gets
      #puts(line)
      if line.include?("overall alignment rate") and x ==2
        alignment << line.split("%")[0]
        x=0
        tmp=temp[0]
        total[tmp]=genome(temp)
        temp=[]
      elsif line.include?("overall alignment rate") and x == 1
        unrrnarate << line.split("%")[0]
        x=2
        tmp=temp[0]
        unpaired[tmp]=unpair(temp)
        temp=[tmp]
      elsif line.include?("overall alignment rate") and x == 0
        pairrrnarate << line.split("%")[0]
        x=1
        tmp=temp[0].chomp
        paired[tmp]=pair(temp)
        temp=[tmp]
      else
        temp << line unless line.include?("SAM") or line.include?("PBS")
      end
    end
   printer(paired,"summary-paired.txt")
   printer(unpaired,"summary-unpaired.txt")
   printer(total,"summary-alignment.txt")
  rescue NoMethodError
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
