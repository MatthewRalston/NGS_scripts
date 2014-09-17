#!/home/mrals/.rvm/rubies/ruby-2.0.0-p247/bin/ruby
=begin
------------------------------------------------------------------------------
--                                                                          --
--                                MATTHEW RALSTON                           --
--                                                                          --
--                        O T H E R _ T R A C K S . R B                     --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Fall 2014                                                               --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed to process foldchange tables from differential     --
-- expression analyses to produce formatted files for input to Circos       --
-- Circos format:                                                           --
-- chr start end value opts...                                              --
-- C 1234 5678 910 svgO1=blah,svgclass=abc                                  --
--                                                                          --
-- SVGs produced by circos now have the ability to add attributes to the svg--
-- image by prefixing svg to optional attributes that you wish to add       --
-- e.g. svgNNNN=abc would result in the NNNN attribute of the svg object    --
-- (heatmap, point, etc.) having a value of abc.                            --
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

DIRECTORY="/home/mrals/Final/circos/data"
GTFPATH="/home/mrals/Final/reference/CAC.gtf"
SIGBUOH="/home/mrals/Final/circos/data/interactions/hierarchical_butanol.raw"
SIGBA="/home/mrals/Final/circos/data/interactions/hierarchical_butyrate.raw"
POSITIONS="/home/mrals/Final/circos/data/positions.txt"
THRESHOLD=0.5
C=3940880
P=192000


################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################

def main
  dict={}
  c=0; p=0
  gtf=gtfread(GTFPATH)
  Dir.glob(DIRECTORY+"/**/*.pre").each do |path|
    dict[path.split('/')[-1].split('.')[0]]=process(path, gtf)
  end
  genes,c,p=significant(gtf)
  dc=(C/c).floor
  dp=(P/p).floor
  Dir.glob(DIRECTORY+"/**/*.pre").each do |path|
    printer(path.split('.')[0]+".fc",dict[path.split('/')[-1].split('.')[0]],gtf,dc,dp,genes,"foldchange")
    printer(path.split('.')[0]+".pval",dict[path.split('/')[-1].split('.')[0]],gtf,dc,dp,genes,"pval")
  end
  bashcommands
end

def gtfread(filepath)
  gtf={}
  File.open(filepath,'r').each do |line|
    line=line.split
    gtf[line[9].gsub(/"/,'').chomp(';')]={"Start"=>line[3],"End"=>line[4],"Strand"=>line[6]}
  end
  gtf
end

def process(filepath,gtf)
  foldchanges={};time=filepath[/\d+/]
  File.open(filepath,'r').each do |line|
    line=line.chomp.split("\t")
    if gtf[line[0]].nil?
      puts filepath
    end
    line[2] == "NA" ? line[2] = 0.9999 : line[2]=line[2].to_f.round(5)
    foldchanges[line[0]]={"foldchange"=>line[1].to_f.round(4),"pval"=>line[2],"time"=>time,"start"=>gtf[line[0]]["Start"],"end"=>gtf[line[0]]["End"],"strand"=>gtf[line[0]]["Strand"],"id"=>line[0]}
  end
  foldchanges
end

def significant(gtf)
  buoh=[];ba=[];genes=[];c=0;p=0
  File.open(SIGBA,'r').each {|line| ba << line.chomp.split[0]}
  File.open(SIGBUOH,'r').each {|line| buoh << line.chomp.split[0]}
  gtf.each {|id,d| genes << id if (ba.include?(id) || buoh.include?(id)) && !id.include?("Ct") && !id.include?("Cr")}
  genes.each {|x| x.include?("P") ? p+=1 : c+=1}
  dc,dp=(C/c).floor,(P/p).floor
  y=0; z=0
  File.open(POSITIONS,'w') do |file|
    genes.each do |x|
      if x.include?("P")
        y+=1
        newline=[x,(y-1)*dp,(y*dp)-1]
      else
        z+=1
        newline=[x,(z-1)*dc,(z*dc)-1]
      end
      file.puts(newline.join("\t"))
    end
  end
  [genes, c, p]
end

def printer(filepath, dictionary, gtf,dc,dp, liszt,value)
  File.open(filepath,'w') do |file|
# x and cd are counters and distances for Chromosomal genes,
# to make the heatmaps equidistant in the svg
    x=0; y=0
    liszt.each do |id|
      value == "pval" ? val=Math.log(dictionary[id][value]/(1-dictionary[id][value])).round(4) : val=dictionary[id][value]
      if id.include?("P")
        y+=1
        newline=["P",(y-1)*dp,(y*dp)-1,val].join(' ')
      else
        x+=1
        newline=["C",(x-1)*dc,(x*dc)-1,val].join(' ')
      end
      liszt=[]
      dictionary[id].keys.each do |key|
        liszt << "svg#{key}=#{dictionary[id][key]}"
      end
      newline+=" #{liszt.join(',')}"
      cond=filepath.split("/")[-1].split(".")[0]
      cond.include?("v") ? newline+=",svgcond=#{cond}" : newline+=",svgcond=#{cond[/\D+/]}"
      file.puts(newline)
    end
  end
end

def bashcommands
  ["fc","pval"].each do |x|
    ["butanol", "butyrate"].each do |y|
      `cat circos/data/foldchanges/stress/#{y}*.#{x} > circos/data/heatmap_#{y}.#{x}`
    end
    `cat circos/data/foldchanges/time/typea/*.#{x} > circos/data/heatmap_time_ti_t1.#{x}`
    `cat circos/data/foldchanges/time/typeb/*.#{x} > circos/data/heatmap_time_ti_ti-1.#{x}`
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
