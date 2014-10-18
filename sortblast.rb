#!/home/mrals/.rvm/rubies/ruby-2.1.2/bin/ruby
=begin
------------------------------------------------------------------------------
--                                                                          --
--                                MATTHEW RALSTON                           --
--                                                                          --
--                           S O R T   B L A S T . R B                      --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Summer 2013                                                             --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed to parse tabular blast results (outfmt6), 
-- separating into groups of exact (one best match output) and multiple     --
-- output.
--                                                                          --
------------------------------------------------------------------------------
=end

################################################
#
#               R E Q U I R E
#
################################################

require 'bio'



################################################
#
#               U S E R    V A R I A B L E S
#
################################################
WORKDIR="Trinity_paired"
#WORKDIR="Trinity_ref"
FASTA="#{WORKDIR}/Trinity-GG.fasta"
BLAST="#{WORKDIR}/Trinity-GG.blast"
DIFF=30
ATTR=[" transcript-length="," align-length="," gap-length="," gaps="," mismatch="," pid="," eval="]
UOUT="#{WORKDIR}/transcripts-unique.gtf"
POUT="#{WORKDIR}/transcripts-partial.gtf"
MOUT="#{WORKDIR}/transcripts-multiple.gtf"


################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################

def parse(line,size)
# transcript  chrom  percent-id  alignment-length  mismatches  gaps  tstart  tend  hitstart  hitend  evalue   score
#                      v
# chrom   source   feature   start   end  score  strand   .   attributes
  line=line.split
  if line[8].to_i > line[9].to_i
    s,e=line[9].to_i,line[8].to_i
    strand="-"
  else
    s,e=line[8].to_i,line[9].to_i
    strand="+"
  end
  chrom,pid,length,mismatch,gaps=line[1..5]
  evalue=line[-2]
  score=line[-1]


  gapsize= e - s - length.to_i + 1
  attributes=[size,length,gapsize,gaps,mismatch,pid,evalue]
  defs=[chrom,"Trinity-blast\ttranscript",s,e,score,strand,".",attributes]
  defs
end

def process(entries,size,unique,multiple,partial)
  name=entries[0].split[0]
  temp=[]
  entries.each {|line| temp << parse(line,size)}
  sorted=temp.sort {|x,y| y[5].to_f<=>x[5].to_f}

  maxes=sorted.select {|x| x[-1]==sorted[0][-1]}
# if there is only one entry in the list (unique match for the alignment of the transcript) add to the unique group
  if maxes.size == 1
    diff=size-(maxes[0][3]-maxes[0][2])
    if diff > DIFF
      partial[name]=sorted
      #puts("#{name} #{diff} #{size} #{maxes[0].to_s}")
      #puts(sorted.to_s)
      #sorted.each {|x| puts(x.to_s)}
    else
      unique[name]=maxes
    end
  else
# if there is more than one entry in the list, we select each element with equal score to the maximum
    multiple[name]=maxes
  end
  [unique,multiple,partial]
end

def iterate(infile,fastas)
# This subroutine gathers blast results from one transcript, passes over all results
# and passes the optimal result(s) to a processing subroutine.
# The results are then added into the unique or multiple dictionaries
  unique={}; multiple={}; partial={}
  temp=[]; transcript=''
  File.open(infile,'r').each do |line|
    #puts(line)
    transcript=line.split[0]
    if temp.any? {|x| x.include?(transcript)}
      temp << line
    else
      if temp.size > 0
        name=temp[0].split[0]
        unique,multiple,partial=process(temp,fastas[name].size,unique,multiple,partial)
      end
      temp = [line]
    end
  end
  puts("unique,partial,multiple")
  puts(unique.size,partial.size,multiple.size)

  [unique, partial, multiple]
end

def gtfout(outfile,dict)
  File.open(outfile,'w') do |file|
    dict.each do |key,value|
      value.each do |liszt|
        line=liszt[0...-1].join("\t")+"\tgeneid=#{key};" + ATTR.zip(liszt[-1]).flatten.join('')
        file.puts(line)
      end
    end
  end
end


def main
  fastas={};Bio::FastaFormat.open(FASTA).each_entry{|f| fastas[f.entry_id]=Bio::Sequence::NA.new(f.seq)}

  unique,partial,multiple=iterate(BLAST,fastas)
  gtfout(UOUT,unique)
  gtfout(POUT,partial)
  gtfout(MOUT,multiple)
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
