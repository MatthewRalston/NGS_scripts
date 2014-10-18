#!/home/mrals/.rvm/rubies/ruby-2.1.2/bin/ruby
=begin
------------------------------------------------------------------------------
--                                                                          --
--                                MATTHEW RALSTON                           --
--                                                                          --
--                               S K E L E T O N . R B                      --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Summer 2013                                                             --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed to adjust a protein fasta and gff3 file, and a     --
-- transcript gtf file to produce the exact coordinates for the protein.    --
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

PROTEINS="Trinotate/Trinity-GG.fasta.transdecoder.pep"
GENOME="reference/CAC.txt"
CDSES="Trinotate/Trin_cdses.gff3"
TRANSCRIPTS="/home/mrals/Final/Trinity_ref/Trinity_filtered.gtf"
R=30

################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################

def gtfread(infile)
  liszt=[]
  File.open(infile,'r').each {|line| liszt << line.split}
  liszt.map! {|x| x[3]=x[3].to_i; x[4]=x[4].to_i; x}
end

def adjustment(trans,peps,protseq,chromseq)
  failures=[];successes=0
  success=[];failtest=[]
  partial=0
  terminators=["taa","tag","tga"]
  transcripts={}; trans.each {|x| transcripts[x[9]]={"chrom"=>x[0],"start"=>x[3],"end"=>x[4],"strand"=>x[6]}}
  peps.each_index do |i|
    record=peps[i]
    unless transcripts[record[0]]
      failures << 'SKIPPING: #{record.join("\t")}'
      next
    end
    #puts(record.join("\t"))
    s=transcripts[record[0]]["start"]
    e=transcripts[record[0]]["end"]
    transcript=chromseq[transcripts[record[0]]["chrom"]][s-1..e-1]
    transcript=transcript.complement if transcripts[record[0]]["strand"] == "-"
    transcript=transcript
    id=record[-1].split(";")[0].split("=")[1][4..-1]
    puts(id)
    protein=protseq[id]
    pstart=record[3]-(R/2)
    pend=record[4]-(R/2)
    if pstart >= transcript.size || pend <= 0
      failures << "#{record}"
      next
    end
    test=successes

    #thing=transcript[record[3]-5..record[3]+5]
    #puts("hooray") if thing.include?("atg")
    #puts thing
    #exit
    temp=[]
    R.times do |x|
      seq=transcript[pstart+x..pend+x]
      #puts(protein.size)
      #puts(seq)
      #puts(seq.translate)
      #puts(protein)
      #exit
      if seq[0..2]=="atg" 
        s=pstart+x
        e=pend+x
        
        if (terminators.include?(seq[-3..-1]) || seq.translate == protein)
          peps[i][3]=peps[i][3]+s
          peps[i][4]=peps[i][4]+e
          puts("HOORAY!!")
          successes+=1
          break
        else
          temp << [s,e]
        end
      end
      failtest << [seq.translate,protein]
    end
    if temp.size > 0
      partial+=1
    end
  end
  puts("SUCCESSES: #{successes}")
  puts("PARTIAL: #{partial}")
  #puts(success)
  puts("FAILURES:")
  puts(failures)
end

def main
  proteins={}; Bio::FastaFormat.open(PROTEINS).each_entry {|f| proteins[f.entry_id]=f.seq}
  fastas={}; Bio::FastaFormat.open(GENOME).each_entry {|f| fastas[f.entry_id]=Bio::Sequence::NA.new(f.seq)}
  transdecoder=gtfread(CDSES)
  transcripts=gtfread(TRANSCRIPTS)
  transcripts.map! {|x| x[9]=x[9].chomp(";").gsub(/"/,''); x}
  adjustment(transcripts,transdecoder,proteins,fastas)
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
