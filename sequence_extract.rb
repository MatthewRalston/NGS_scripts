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
-- This file is designed to extract genomic sequences (nt) from a genome    --
-- given a gtf as input.
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
GTF,FASTA,OLD=ARGV[0..2]
usage='sequence-extract.rb annotation.gtf genomic-sequence.fasta'
puts(usage) if GTF=="-h"




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
    gtf << line.chomp.split
  end
  gtf.map! {|x| x[3]=x[3].to_i;x[4]=x[4].to_i;x[9]=x[9].chomp(";").gsub(/"/,'');x}
end

def extract(gtf,fastas)
  seqs={}
  chrom,plas=fastas.keys
  gtf.each do |gene|
    seq=(gene[0] == chrom ? fastas[chrom][gene[3]-1..gene[4]-1] : fastas[plas][gene[3]-1..gene[4]-1])
    seq=(gene[6] == "+" ? seq : seq.complement)
    seqs[gene[9]]=seq
  end
  seqs
end


def fastaout(fastas)
  fastas.each {|id,seq| STDOUT.puts(">#{id}\n#{seq.upcase}")}
end

def main
  gtf=gtfread(GTF)
  fastas={}; Bio::FastaFormat.open(FASTA).each_entry {|f| fastas[f.entry_id]=Bio::Sequence::NA.new(f.seq)}
  seqs=extract(gtf,fastas)
  oldfastas={}; Bio::FastaFormat.open(OLD).each_entry {|f| oldfastas[f.entry_id]=Bio::Sequence::NA.new(f.seq)}
  test=seqs.keys[0]
  #puts("#{test}:\n#{seqs[test]}\nvs:\n#{oldfastas[test]}")
  fastaout(seqs)
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
