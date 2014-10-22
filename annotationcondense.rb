#!/home/mrals/.rvm/rubies/ruby-2.1.2/bin/ruby
=begin
------------------------------------------------------------------------------
--                                                                          --
--                                MATTHEW RALSTON                           --
--                                                                          --
--               A N N O T A T I O N - C O N D E N S E . R B                --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Summer 2013                                                             --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed to condense three annotations together, conserving --
-- the maximal amount of annotation information. Annotations with identical --
-- coordinates are merged into a single consensus entry.
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

WORKDIR="Trinity_paired"
# Standard cds was preprocessed to join
#cat Trinity_paired/standard-cds.gtf | ruby -ne 'line=$_.split("\t"); line[-1]=line[-1].split(";").map{|x| x.split.join("=")}.join("; "); puts(line.join("\t"))' > Trinity_paired/standard-cds.equals.gtf
STANDARD="#{WORKDIR}/standard-cds.equals.gtf"
TRANSDECODER="#{WORKDIR}/transdecoder.rast.gtf"
RAST="#{WORKDIR}/rast.output.gtf"
OUTPUT="#{WORKDIR}/consensus.cds.gtf"
Chrom="NC_003030.1"

################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################

def gtfread(infile)
  gtf=[]
  File.open(infile,'r').each do |line|
    if line[0] != "#"
      line= line.split("\t")
      line[-1]=line[-1].chomp.split(";").map{|x|x.split(" ").join("_")+";"}
      gtf << line.flatten
    end
  end
  gtf.map!{|x| x[3]=x[3].to_i; x[4]=x[4].to_i; x}
end

def merge(liszt)
  if liszt.size == 1
    return liszt
  else
    id=[]; locus=[]; db=[]; parent=[]; name=[]; ids=[]
    liszt.each do |record|
      record.each do |x|
        y=x.to_s
        id << y if y.include?("ID=")# && y.include?("CA_")
        ids << y if y.include?("_id")
        name << y if y.include?("Name")

        locus << y if y.include?("locus_tag")
        db << y if y.include?("db_xref")
        parent << y if y.include?("Parent")
      end
    end
  end
  i=id.index{|x| x.include?("CA_")}
  id=[(i ? id[i] : id[0])]
  i=liszt.index{|x|x[1]=="gb2gtf"}
  record=(i ? liszt[i] : liszt[0])
  name=(name.size > 1 && name[0]!=name[1]) ? ["Name="+name.map{|x|x.split("=")[1].chomp(";")}.join(",")+";"] : [name[0]]
  extras=(id+locus+parent+ids+db+name).select{|x|x}
  extras.map!{|x| x[-1]==";" ? x : x+";" }
  liszt=record[0..7].join("\t")+"\t"
  liszt+=" "+extras.join(" ")
  liszt=liszt.split
  liszt[3]=liszt[3].to_i; liszt[4]=liszt[4].to_i
  liszt
end

def condense(ref,gtf)
  # create a combined list to iterate through
  combined=(ref+gtf)
  combined=(combined.select{|x|x[0]==Chrom}.sort{|x,y|x[3]<=>y[3]}+combined.select{|x|x[0]!=Chrom}.sort{|x,y|x[3]<=>y[3]})
  # create empty list for related records and condensed records
  group=[]; condensed=[]; test=[0,3,4,6]; l=[0,3,6]; r=[0,4,6]
  i=0
  while record=combined[i]
    if group.size > 0 && test.map{|x|record[x]} == test.map{|x|group[0][x]}
      group << record
      group=[merge(group)]
    elsif group.size > 0 && (l.map{|x|record[x]}==l.map{|x|group[0][x]} || r.map{|x|record[x]}==r.map{|x|group[0][x]})
      group << record
      s=group.map{|x| x[3]}.sort[0]
      e=group.map{|x|x[4]}.sort[-1]
      group.map!{|x|x[3]=s;x[4]=e;x}
      group=[merge(group)]
    else
      condensed << group[0] if group.size > 0
      raise "multiple entries #{group.to_s}" if group.size > 1
      group = [record]
    end
    i+=1
  end
  puts(condensed.size)
  condensed
end

def gtfout(gtf,outfile)
  File.open(outfile,'w') do |file|
    gtf.each {|x| file.puts(x[0..7].join("\t")+"\t"+x[8..-1].join(" "))}
  end
end

def main
  standard=gtfread(STANDARD)
  transdecoder=gtfread(TRANSDECODER)
  rast=gtfread(RAST)
  puts("standard:#{standard.size}\ntransdecoder:#{transdecoder.size}\nrast:#{rast.size}")
  final=condense(standard,transdecoder)
  final=condense(final,rast)
  gtfout(final.sort{|x,y|x[3].to_i<=>x[4].to_i},OUTPUT)
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
