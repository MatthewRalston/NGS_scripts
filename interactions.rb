#!/home/mrals/.rvm/rubies/ruby-2.0.0-p247/bin/ruby
=begin
------------------------------------------------------------------------------
--                                                                          --
--                                MATTHEW RALSTON                           --
--                                                                          --
--                           I N T E R A C T I O N . R B                    --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Summer 2013                                                             --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed to process interaction records into circos input   --
-- format.
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
GTFPATH="/home/mrals/Final/reference/CAC.gtf"
CLUSTERS="/home/mrals/Final/circos/data/interactions/hierarchical_butyrate.raw"
OUTFILE="/home/mrals/Final/circos/data/interactions/hierarchical_butyrate.txt"
#CLUSTERS="/home/mrals/Final/circos/data/interactions/hierarchical_butanol.raw"
#OUTFILE="/home/mrals/Final/circos/data/interactions/hierarchical_butanol.txt"
LABELS="/home/mrals/Final/circos/data/genelabels.txt"
POSITIONS="/home/mrals/Final/circos/data/positions.txt"
# 10 color scheme - Diverging, half up
#COLORS={"1"=>"(165,0,38)","2"=>"(215,48,39)","3"=>"(244,109,67)","4"=>"(253,174,97)","5"=>"(254,224,139)",
#  "6"=>"(217,239,139)","7"=>"(166,217,106)","8"=>"(102,189,99)","9"=>"(26,152,80)","10"=>"(0,104,55)"}
# 8 color scheme - Diverging
#COLORS={"1"=>"(165,0,38)","2"=>"(215,48,39)","3"=>"(244,109,67)","4"=>"(253,174,97)",
#  "5"=>"(166,217,106)","6"=>"(102,189,99)","7"=>"(26,152,80)","8"=>"(0,104,55)"}
# 6 color scheme - Diverging
#COLORS={"1"=>"(165,0,38)","2"=>"(215,48,39)","3"=>"(244,109,67)",
#  "4"=>"(102,189,99)","5"=>"(26,152,80)","6"=>"(0,104,55)"}
# 4 color scheme - Diverging
#COLORS={"1"=>"(165,0,38)","2"=>"(215,48,39)",
#  "3"=>"(26,152,80)","4"=>"(0,104,55)"}
# 11 color scheme - Clasification

COLORS={
"1"=>"one",
"2"=>"two",
"3"=>"three",
"4"=>"four",
"5"=>"five",
"6"=>"six",
"7"=>"seven",
"8"=>"eight",
"9"=>"nine",
"10"=>"ten",
"11"=>"eleven",
"12"=>"twelve",
"13"=>"thirteen",
"14"=>"fourteen",
"15"=>"fifteen",
"16"=>"sixteen",
"17"=>"seventeen",
"18"=>"eighteen",
"19"=>"nineteen",
"20"=>"twenty",
"21"=>"twentyone",
"22"=>"twentytwo",
"23"=>"twentythree",
"24"=>"twentyfour"}

N=COLORS.size
################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################

def main
  c=[]
  File.open(CLUSTERS,'r').each {|line| c << line.chomp.split("\t") unless (line.include?("Cr") || line.include?("Ct"))}
  gtf=positions
  process(c,gtf)
  printlabels(gtf)
end

def positions
  gtf={}
  File.open(POSITIONS,'r').each do |line|
    line=line.split
    gtf[line[0]]={"Start"=>line[1],"End"=>line[2]}
  end
  gtf
end

def process(liszt,gtf)
  temp={}; first=''
  # create a temporary storage for each cluster with genes in a list and a count of the
  # number of these genes
  (1..N).to_a.each {|x| temp[x.to_s] = {"chr"=>[],"plas"=>[],"p"=>0,"c"=>0}}
  # add genes to the lists, incrementing the counter
  liszt.each {|gene,c| (gene.include?("P") ? (temp[c]["plas"] << gene;temp[c]["p"]+=1) : (temp[c]["chr"] << gene; temp[c]["c"]+=1))}
  tmp=[]
  temp.each {|id,d| tmp << id if !tmp.include?(id) && (d["c"] == 1)}
  puts tmp
  File.open(OUTFILE,'w') do |file|
    temp.each do |cluster,d|
      firstp,firstc=d["plas"][0],d["chr"][0]
      if d["c"] > 2
        (d["c"]-1).times do |x|
          gene=d["chr"].shift
          file.puts(link("C",gene,d["chr"][0],gtf,cluster))
        end
        file.puts(link("C",d["chr"][0],firstc,gtf,cluster))
      elsif d["c"] == 2
        file.puts(link("C",firstc,d["chr"][1],gtf,cluster))
      end
      if d["p"] > 2
        (d["p"]-1).times do |x|
          gene=d["plas"].shift
          file.puts(link("P",gene,d["plas"][0],gtf,cluster))
        end
        file.puts(link("P",d["plas"][0],firstp,gtf,cluster))
      elsif d["p"] == 2
        file.puts(link("P",firstp,d["plas"][1],gtf,cluster))
      end
    end
  end
end

def link(c,g1,g2,gtf,i)
  [c,gtf[g1]["Start"],gtf[g1]["End"],c,gtf[g2]["Start"],gtf[g2]["End"],"svgid=#{g1},svgpartner=#{g2},svgclass=Cluster-#{i},color=#{COLORS[i]}"].join(" ")
end

def printlabels(gtf)
  File.open(LABELS,'w') do |file|
    gtf.each do |gene,d|
      if gene.include?("P")
        newline=["P",d["Start"],d["End"],gene,"svgid=#{gene}"]
      else
        newline=["C",d["Start"],d["End"],gene,"svgid=#{gene}"]
      end
      file.puts(newline.join(" "))
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
