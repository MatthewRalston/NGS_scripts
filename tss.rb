#!/home/mrals/.rvm/rubies/ruby-2.0.0-p247/bin/ruby
=begin
------------------------------------------------------------------------------
--                                                                          --
--                                MATTHEW RALSTON                           --
--                                                                          --
--                                   T S S . R B                            --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Summer 2013                                                             --
--                                                                          --
------------------------------------------------------------------------------
-- This script is designed to locate and summarize transcription start sites--
-- 1. A gtf file is loaded into memory and a window near the start of each  --
--      gene is written to a BED format file.                               --
-- 2. BEDtools BAMtoBED is used to convert the file to BED format, piped to --
-- 3. The BED file is then sorted and passed to BEDtools intersect          --
-- 4. BEDtools intersect is used to find the set intersections of the       --
--      reads with the transcription initiation region, written to BED files--
-- 5. The BED file will then be sorted according the coordinates of the     --
--      reads and then read into memory.                                    --
-- 6. For each tss, its corresponding records will be loaded into a list.   --
-- 7. This list of BED records is then passed to BEDtools genomecov.        --
-- 8. The resulting bedgraph output is returned as a list and processed     --
-- 9. The processed bedgraph is suitable for processing by TSSi using       --
--      the gem 'rinruby.' After proceeding through regularization and      --
--      identification, a two item vector is returned containing the TSS    --
--      and the normalized expression value at this location.               --
-- 10. Create plot of TSS according to cufflinks
--                                                                          --
------------------------------------------------------------------------------
=end

################################################
#
#               R E Q U I R E
#
################################################

require 'rserve'


################################################
#
#               U S E R    V A R I A B L E S
#
################################################

CUFFBED="Cuffcompare/canonical.combined.bed"
# Quick pipe to convert a gtf into BED6
# cat myfile | gtf2bed | bed12ToBed6 > mynewfile
# to filter add ->    | egrep -v "exclude" |
REFGENOME="/home/mrals/ETP/reference/CAC.genome"
# This is a list of strings of gene_ids (e.g. CA_C2020) to plot with ggplot2
TOPLOT=["CA_P0001"]
PLOTDIR="/home/mrals/ETP/tsses/"
TIRFILE="/home/mrals/ETP/tsses/tir.bed"
RADIUS=300
CON=Rserve::Connection.new
CON.eval("library('TSSi')")
CON.eval("library('ggplot2')")
CON.eval('setwd("/home/mrals/ETP/")')
BAMS=["SAM_processed/NS30A.3.bam"]
#BAMS=`ls SAM_processed/*.3.bam`.split("\n")

################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################

def ref_read(infile)
# This function reads the reference annotation bed file into a list
  liszt=[]
  File.open(infile,'r').each do |line|
    liszt << line.split("\t")
  end
  liszt
end

def ref_mod(annotation,radius,outfile)
# This function reads a 2D list of BED format annotations and prints a modified BED file of the TIRs around the
# annotation 'start' site with a window defined by [start-radius : start+radius]
# This function returns this 2D list of TIRs in BED format.
  liszt=[]
  File.open(outfile,'w') do |file|
    annotation.each do |l|
      l[5] == "+" ? temp = [l[0],l[1].to_i-radius,l[1].to_i+radius]+l[3..8] : temp = [l[0],l[2].to_i-radius,l[2].to_i+radius]+l[3..8] 
      temp[1]=0 if temp[1] < 0
      liszt << temp
      file.puts((temp).join("\t"))
    end
  end
  liszt
end

def intersection(bamfile, tirfile, annotation)
# THE BIG ONE
# This function takes the name of a BAM file, a gtf format TIR file, and a BED format list of gene annotation.
# This function returns a hash (by gene_id) of tsses and their expression levels.
# This function intersects a bamfile with a BED file of transcription initiation regions of each gene
# Then for each gene in the 'annotation' list, the reads from that area are found and passed to tssfinder
  x=0
  tss={}
  reads=[]
  tempfile="tsses/temp"
  `samtools view -bh #{bamfile} | bam2bed | ruby -ne 'puts($_.split("\t")[0...6].join("\t"))' | intersectBed -a stdin -b #{tirfile} | sort -k 1,1 -k2,2n > #{tempfile}`
  File.open(tempfile,'r').each {|read| reads << read.chomp.split("\t")}
  annotation.each do |liszt|
# here each record of the annotation is used to locate the reads in this are with binary search, and pass them
# and the annotation to the function tssfinder
    a,b=(0...reads.size).to_a.bsearch{|c| reads[c][1].to_i >= liszt[1].to_i},(0...reads.size).to_a.bsearch{|c| reads[c][2].to_i > liszt[2].to_i}
    unless a==nil or b==nil
      t1,t2=tssfinder(reads[a...b],liszt,bamfile,x)
      tss[liszt[3]]={'tss'=>t1,'expr'=>t2}
    end
    x+=1
  end
  tss
end

def tssfinder(reads,annotation,bam,x)
# ON EACH GENE
# This function is passed a 2D array of BED format reads, split along tabs
# It is also given a list of a BED format record.
# First, coverage is calculated for this list, then the results are processed and written to file.
# Next, a function is executed in R that reads this list and statistically evaluates the transcription start site
  tempfile="tsses/tmp"
  bg=[]
  File.open(tempfile,'w') {|file| reads.each {|read| file.puts(read.join("\t"))}}
  `bedtools genomecov -bg -i #{tempfile} -g #{REFGENOME}`.split("\n").each do |bed|
    bed=bed.split("\t")
    bg << [bed[0], 1, bed[1], annotation[5], bed[-1]]
  end
  File.open(tempfile,'w') {|file| bg.each {|record| file.puts(record.join("\t")) if (record[2].to_i >= annotation[1].to_i and record[2].to_i <= annotation[2].to_i)}}
  tsses,reads=tssi(tempfile,annotation,bam,x)
end

def tssi(infile,annotation,bam,x)
# This function runs a tss prediction for each file 
  input=CON.eval("input<-read.table('#{infile}',header=F)")
  CON.eval("attach(input)")
  CON.eval("x<-segmentizeCounts(chr=V1,region=V2,start=V3,strand=V4,counts=V5)")
  CON.eval("detach(input)")
  CON.eval("my.ratio<-normalizeCounts(x)")
  CON.eval("my.fit<-normalizeCounts(x,fit=TRUE)")
  CON.eval("res<-identifyStartSites(my.fit)")
  # A ggplot plot will produced of local coverage and the computed TSS if the id is present in the
  # TOPLOT list
  tssplot(annotation,bam,infile,x) if TOPLOT.include?(annotation[3])
  tsses,reads=CON.eval("tsses<-tss(res)").to_ruby.to_a
end

def tssplot(annotation,bam,tempfile,x)
# This function takes a BED format annotation list, the name of a BAM format file, and the name of a temporary file
# This function returns nothing.
# This function produces 
  id=annotation[3]
  covfinder(bam,annotation,tempfile)
  CON.eval('t<-tss(res,1)')
  df=CON.eval("df<-reads(res,1)")
  p=CON.eval("p<-ggplot()+geom_area(aes(x=cov$V2,y=cov$V3,fill='Reads'))+geom_area(data=df,aes(x=start,y=fit,fill='Fit'))+geom_segment(aes(x=t$pos,y=rep(0,length(t$pos)),xend=t$pos,yend=t$reads))")
  CON.eval("ggsave('#{PLOTDIR}#{id}_#{x}.png',p)")
end

def covfinder(bamfile,annotation,tempfile)
# This function calculates local coverage for each base in the annotation's window
# Next, this function reads in the list to a dataframe object in R
  `samtools view -bh #{bamfile} "#{annotation[0]}:#{annotation[1].to_i+1}-#{annotation[2].to_i+1}" | bam2bed | ruby -ne 'puts($_.split("\t")[0...6].join("\t"))' | bedtools genomecov -d -strand #{annotation[5]} -i stdin -g #{REFGENOME} | ruby -ne 'puts $_ if ($_.split[1].to_i > #{annotation[1].to_i+1} and $_.split[1].to_i <= #{annotation[2].to_i+1})' > #{tempfile}`
  CON.eval("cov<-read.table('#{tempfile}',header=F)")
end

def main
  raw=ref_read(CUFFBED)
  tir=ref_mod(raw,RADIUS,TIRFILE)
  BAMS.each do |bamfile|
    tss=intersection(bamfile,TIRFILE,tir)
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
