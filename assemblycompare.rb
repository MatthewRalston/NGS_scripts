#!/home/mrals/.rvm/rubies/ruby-2.0.0-p247/bin/ruby
=begin
------------------------------------------------------------------------------
--                                                                          --
--                                MATTHEW RALSTON                           --
--                                                                          --
--                A S S E M B L Y   C O M P A R E. R B                      --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Summer 2013                                                             --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed to identify items in a first annotation that are   --
-- completely contained in genetic elements of another.                     --
-- other items of the list. This is meant to find the transcripts that are  --
-- fully contained in other transcripts.                                    --
-- Additionally, this script produces some statistics about the assembly.   --
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



GTF,REFERENCE=ARGV[0..1]
OUTDIR=ARGV[2]
NOVEL="#{OUTDIR}/novel.gtf"
STANDARD="#{OUTDIR}/standard.gtf"
OLD="#{OUTDIR}/standard-cds.gtf"
OLDGENESLIST="#{OUTDIR}/standard-cds-list.txt"
STATS="#{OUTDIR}/stats"
DIST=40

################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################


def gtfread(infile)
  liszt=[]
  File.open(infile,'r').each {|line| liszt << line.chomp.split}
  liszt.map! {|x| x[3]=x[3].to_i; x[4]=x[4].to_i; x}
  liszt.sort! {|x,y| x[3]<=>y[3]}
end

def contains?(l1,l2)
  # if second contains first
  test=(l2[3] <= l1[3] && l1[4] <= l2[4] && l1[0] == l2[0] && l1[6] == l2[6])
  test
end


def contain(l1,l2)
  # if reference(L1) is inside new annotation (L2), add L2 to 
  contains=[]; contained=[]
  l1.each do |g1|
    l2.each do |g2|
      if contains?(g1,g2)
        contains << g2
        contained << g1
        break
      end
    end
  end
  [contains,contained]
end

def gtfout(records,output)
  File.open(output,'w') do |file|
    records.each do |x|
      file.puts(x.join("\t"))
    end
  end
end

def geneids(records,output)
  File.open(output,'w') do |file|
    records.each do |x| 
      if x[12]     
        file.puts(x[-1].split(":")[1].gsub(/"/,""))
        next
      else
        file.puts(x[9].chomp(";").gsub(/"/,""))
        next
      end
    end
  end
end

def tlength_distr(gtf,outfile)
  File.open(outfile,'w') do |file|
    gtf.each {|rec| file.puts(rec[4]-rec[3])}
  end
end

def listout(liszt,outfile)
  File.open(outfile,'w') {|file| liszt.each {|x| file.puts(x)}}
end

def operon(transcripts,cds)
  dictionary={}
  tlength=[]
  olength=[]
  five=[]
  three=[]
  transcripts.each do |t|
    liszt=[]
    cds.each {|c| liszt << c if contains?(c,t)}
    if liszt.size > 0
      liszt.sort! {|x,y| x[3]<=>y[3]}
      # If + strand, min CDS start - t[start]   :  t[end] - max CDS end
      five << (t[6] == "+" ? liszt.map{|x| x[3]}.min-t[3] : t[4]-liszt.map{|x|x[4]}.max)
      # If + strand, t[end] - max CDS end     : min CDS start - t[start]
      three << (t[6] == "+" ? t[4]-liszt.map{|x|x[4]}.max : liszt.map{|x|x[3]}.min-t[3])
      if liszt.size > 1
        tlength << t[4]-t[3]
        olength << liszt.size
      end
    end
    dictionary[t[9]]=liszt
  end
  listout(tlength,STATS+"/operon.len")
  listout(olength,STATS+"/operonsize.txt")
  listout(five,STATS+"/five.len")
  listout(three,STATS+"/three.len")
end

def igrlength_distr(liszt,outfile)
  igrs=[]
  old=nil
  liszt.each do |plus|
    if plus[6] == "+"
      if old && plus && old[6] == "+"
        igrs << plus[3]-old[4] if plus[3] > old[4]
      end
      old=plus
    end
  end
  old=nil
  liszt.each do |minus|
    if minus[6] == "-"
      if old && minus && old[6] == "-" && minus[6] == "-"
        igrs << minus[3]-old[4] if minus[3] > old[4]
      end
      old=minus
    end
  end
  listout(igrs,outfile)
end

def antisense(liszt,outfile)
  pairs=[]
  liszt.each do |t1|
    liszt.each do |t2|
      left=(t2[4] > t1[3] && t2[4] < t1[4])
      right=(t2[3] > t1[3] && t2[3] < t1[4])
      if t1[6] != t2[6] && (left || right)
        left=[t1[3],t2[3]].max
        right=[t1[4],t2[4]].min
        l1,l2=[t1[4]-t1[3],t2[4]-t2[3]]
        if (l1 < 400 || l2 < 400) && right-left >= DIST
          one,two=[t1[9].chomp(";").gsub(/"/,''),t2[9].chomp(";").gsub(/"/,'')]
          new=[]
          l1 > l2 ? new=[two,one] : new=[one,two]
          pairs << new unless pairs.include?(new)
        end
      end
    end
  end
  File.open(outfile,'w') {|file| pairs.each {|x| file.puts(x.join("\t"))}}
end

def stats(total,novel,reference,cds,refcds)
  `mkdir #{OUTDIR}/stats`
  tlength_distr(novel,STATS+"/novel.len")
  tlength_distr(reference,STATS+"/standard.len")
  tlength_distr(refcds,STATS+"/refcds.len")
  operon(reference,cds)
  igrlength_distr(total,STATS+"/igr.len")
  igrlength_distr(refcds,STATS+"/refigr.len")
  antisense(total,STATS+"/antisense.pairs")
end

def main
  gtf=gtfread(GTF)
  reference=gtfread(REFERENCE)
  contains,contained=contain(reference,gtf)
  novel=gtf.dup.delete_if {|x| contains.include?(x)}
  gtfout(novel,NOVEL)
  gtfout(gtf-novel,STANDARD)
  gtfout(contained,OLD)
  #geneids(contained,OLDGENESLIST)
  stats(gtf,novel,gtf-novel,contained,reference)
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
