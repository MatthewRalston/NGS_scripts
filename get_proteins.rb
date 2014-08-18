#!/home/mrals/.rvm/rubies/ruby-2.0.0-p247/bin/ruby
# This script acquires proteins from a gtf, printing the proteins to STDOUT in fasta format

require 'bio'

REFGENOME="reference/CAC.txt"

def fasta_read(infile)
  fastas={}
  Bio::FastaFormat.open(infile).each_entry {|f| fastas[f.entry_id]=Bio::Sequence::NA.new(f.seq)}
  fastas
end



def gene_receive
  genes={}
  while line=STDIN.gets
    line=line.split
    genes[line[9][1...-2]]={"chr"=>line[0],"start"=>line[3].to_i,"end"=>line[4].to_i,"strand"=>line[6]}
  end
  genes
end

def prot_print(genes,fastas)
  genes.each do |key,entry|
    if entry["strand"] == "+"
      protein=fastas[entry["chr"]][entry["start"]-1...entry["end"]].translate[0...-1]
    else
      protein=fastas[entry["chr"]][entry["start"]-1...entry["end"]].complement.translate[0...-1]
    end
    STDOUT.puts(">#{key}\n#{protein}")
  end
end


def main
  genes=gene_receive
  fastas=fasta_read(REFGENOME)
  prot_print(genes,fastas)
end


main
