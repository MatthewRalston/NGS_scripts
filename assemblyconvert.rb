#!/home/mrals/.rvm/rubies/ruby-2.0.0-p247/bin/ruby
=begin
------------------------------------------------------------------------------
--                                                                          --
--                                MATTHEW RALSTON                           --
--                                                                          --
--                 A S S E M B L Y  C O N V E R T . R B                     --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Summer 2013                                                             --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed to convert a gtf annotation of CAC (produced from  --
-- the assembly project) to a genbank file, suitable for RAST.              --
-- Step 1. Load the reference genome fasta, transdecoder gtf, and transcript--
--   alignment gtf.
-- Step 2. convert the transdecoder gff3 into genomic coordinates
-- Step 3. Print the adjusted gtf and run seqret to generate genbank format --
-- Step 4. Load the peptide sequences from the genomic coordinates
-- Step 5. Parse genbank and add the CDS sequences




-- #NOTE: The gtf file must consist of 'gene' entries only and also be an   --
-- appropriately formated gtf file, with the goofy fields at the end        --
-- formatted with spaces. See below:                                        --
cat Trinity_ref/Trinity_filtered.gtf | sed 's/transcript\t/gene\t/' | ruby -ne 'line=$_.chomp.split; puts(line[0..7].join("\t")+"\t"+line[8..-1].join(" "))' | sort -k 1,1 -k4,4n > Trinity_ref/Trin_as_genes.gtf
--
--
--
-- Then, both genbank files must be concatenated into one file
--
--
-- Remember to change the Article title 'ART-TITLE', Journal title          --
-- 'JOURNAL-TITLE', Pubmed id 'PUBMED-ID', and date of production 'X/X/XXXX'--
--
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
ANNOTATEDIR="Trinotate2"
FASTA="reference/CAC.txt"
TRANSDECODER="#{ANNOTATEDIR}/Trin-blast.complete.cds"
TRANSCRIPTS="#{WORKDIR}/Trin-blast.gtf"
TRANSCRIPTFASTA="#{WORKDIR}/Trin-blast.fasta"
CDSFASTA="#{ANNOTATEDIR}/Trin-blast.fasta.transdecoder.complete.cds.fasta"
STOPS=["TAA","TAG","TGA"]
FINALCDS="#{ANNOTATEDIR}/Trin-blast.annotation.gtf"
FINALCHROM="#{ANNOTATEDIR}/Trin-blast.annotation.chrom.gtf"
FINALPLAS="#{ANNOTATEDIR}/Trin-blast.annotation.plas.gtf"
#STANDARDCDS="#{WORKDIR}/standard-cds.gtf"
STANDARDCDS="#{WORKDIR}/consensus.cds.space.gtf"
PREFIX=" "*21
TEMP="combined.genbank"
CHROMFASTA="reference/chrom.ffn"
PLASFASTA="reference/plas.ffn"
CHROMGB="reference/chromgenbank.txt"
PLASGB="reference/plasgenbank.txt"
CHROMFEATURE='      source          1..3940880
                      /organism="Clostridium acetobutylicum ATCC 824"
                      /mol_type="genomic DNA"
                      /strain="ATCC 824"
                      /db_xref="ATCC:824"
                      /db_xref="taxon:272562"'
PLASFEATURE='     source          1..192000
                     /organism="Clostridium acetobutylicum ATCC 824"
                     /mol_type="genomic DNA"
                     /strain="ATCC 824"
                     /db_xref="ATCC:824"
                     /db_xref="taxon:272562"
                     /plasmid="pSOL1"
                     /note="studies by P. Soucaille and coworkers (INSA,
                     Toulose) have shown that loss of this plasmid coincides
                     with the loss of the capacity to produce acetone and
                     butanol and that the genes involved in solvent formation
                     reside on pSOL1"'
TRANSDECODEROUT="#{ANNOTATEDIR}/transdecoder.rast"
MERGEDOUT="#{ANNOTATEDIR}/merged.rast"
STANDARDOUT="#{ANNOTATEDIR}/standard.rast"
Chrom="NC_003030"

################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################

def gtfread(infile)
  liszt=[]
  File.open(infile,'r').each do |line|
    line=line.split
    liszt << line if line[0]
  end
  liszt.map! {|x| x[3]=x[3].to_i; x[4]=x[4].to_i; x}
end

def coordfix(transdecoder,cds,fastas)
  transdecoder.map! do |x|
    l=x[3]
    (l-3 < 0 ? l=0 : l=l-3)
    c=cds[x[-1]]
    seq=fastas[x[-1].split("|")[0]]
    i,j=nil,nil
    3.times{|y| i=y if seq[l+y..l+y+2] == "ATG"}
    3.times{|y| j=y if STOPS.include?(seq[l+c.size+y-3..l+c.size+y-1])}
    raise x.to_s unless i == j
    l,r=l+i,l+c.size+j-1
    x[3]=l
    x[4]=r
    x << seq[r-2..r]
    x
  end
  transdecoder
end

def converter(transcriptlist,transdecoder,genome)
  transcript={}
  transcriptlist.each {|x| transcript[x[9]]={"start"=>x[3],"end"=>x[4],"strand"=>x[6],"chrom"=>x[0]}}
  transdecoder.map! do |feature|
    feature[6]=transcript[feature[0]]["strand"]
    chrom=transcript[feature[0]]["chrom"]
    l=transcript[feature[0]]["start"]-1
    r=transcript[feature[0]]["end"]-1
    #puts("T:#{l},#{r}")
    #puts("CDS:#{feature[3]},#{feature[4]}")
    #puts(feature.to_s)
    seq=genome[chrom][l..r]
    seq=feature[6] == "+" ? seq[feature[3]..feature[4]] : seq[seq.size-feature[4]-1..seq.size-feature[3]-1].complement
    #puts(seq)
    #puts("T in CDS:#{l+feature[3]} - #{l + feature[4]}")
    if feature[6] == "+" && seq[0..2]=="atg" && seq[-3..-1].upcase == feature[-1]
      r=l+feature[4]+1
      l=l+feature[3]+1
    elsif feature[6] == "-" && seq[0..2] == "atg" && seq[-3..-1].upcase == feature[-1]
      l=r-feature[4]+1
      r=r-feature[3]+1
    else
      raise "ERROR: #{feature.to_s}"
    end
    s=genome[chrom][l-1..r-1]
    s=s.complement if feature[6] == "-"
    raise "error: #{feature.to_s}" unless seq==s
    feature[3..4]=[l,r]

    feature.include?("Parent") ? feature[feature.index("Parent")+1]=feature[0] : (i=feature.index("ID"); feature[i+1]=feature[i+1].split("|")[0])
    feature[0]=chrom
    feature[0..-2]
  end
  [transdecoder,transcript]
end

def gtfout(gtf,outfile)
  File.open(outfile,'w') do |file|
    gtf.each do |x|
      line=x[0..7].join("\t")+"\t"+x[8..-1].join(" ")
      line=line.split(" Parent ").join("; Parent=").split("ID ").join("ID=")+";"
      file.puts(line)
    end
  end
end

def seqret(fasta,gff)
 `seqret -sequence #{fasta} -feature -fformat gff -fopenfile #{gff} -osformat genbank -auto`
end

def proteinextract(gtf,fastas)
  peptides={}
  gtf.each do |x|
    seq=fastas[x[0]][x[3]-1..x[4]-1]
    seq=seq.complement if gtf[6] == "-"
    raise "Transcript size error\n#{seq}\nsize#{seq.size}\n#{seq.size%3}" unless seq.size%3==0
    key=x[9]#[4..-1].gsub(/\|/,"-").gsub(/\./,"-")
    peptides[key]=seq.translate
  end
  peptides
end


def readnwrite(input,output,gtf,proteins)
  currentchrom=''
  features=false
  currentfeature=''
  chromrecords=gtf.select{|x| x[0]==Chrom}
  plasrecords=gtf.select{|x| x[0]!=Chrom}
  File.open(output,'w') do |file|
    File.open(input,'r').each do |line|
      if line.include?('LOCUS')
        currentfeature=''
        features=false
        currentchrom = locus(line)
        file.puts(line)
      elsif line.include?('DEFINITION')
        file.puts(line)
        addheader(currentchrom,file)
      elsif line.include?('FEATURES')
        file.puts(line)
        featureprefix(currentchrom,file)
        features=true
      elsif line.include?("CDS") 
        currentfeature="CDS"
        file.puts(line)
      elsif line.include?("gene")
        currentfeature="gene"
        file.puts(line)
      elsif line.include?('cds') && line.include?("locus_tag")
        file.puts(line)
        temp=PREFIX+"/translation=\""
        key=line.chomp.split("=")[1].gsub(/"/,'')
        temp+=proteins[key]+"\""
        file.puts(temp)
      else
        file.puts(line)
      end
    end
  end
end


def locus(chromline)
  chromline.split[1]
end

def addheader(current,file)
  if current.include?("NC_003030")
    file.puts(File.open(CHROMGB,'r').read)
  elsif current.include?("NC_001988")
    file.puts(File.open(PLASGB,'r').read)
  else
    raise "Unknown chromosome in addheader: #{current}"
  end
end

def featureprefix(current,file)
  if current.include?("NC_003030")
    file.puts(CHROMFEATURE)
  elsif current.include?("NC_001988")
    file.puts(PLASFEATURE)
  else
    raise "Unknown chromosome in feature prefix: #{current}"
  end
end

def main
# Step 1. Load the reference genome fasta, transdecoder gtf, and transcript alignment gtf.

  transdecoder=gtfread(TRANSDECODER).map {|x| x[-1]=x[-1].split(/[;=]/);x.flatten}
  transcripts=gtfread(TRANSCRIPTS).map {|x| x[8..9]=x[8].chomp(";").gsub(/"/,'').split("=");x}
  genome={}; Bio::FastaFormat.open(FASTA).each_entry {|f| genome[f.entry_id]=Bio::Sequence::NA.new(f.seq)}
  cds={};Bio::FastaFormat.open(CDSFASTA).each_entry {|f| cds[f.entry_id.gsub(/~/,"|")]=f.seq}
  transeq={};Bio::FastaFormat.open(TRANSCRIPTFASTA).each_entry {|f| transeq[f.entry_id]=f.seq}
# Step 2. convert the transdecoder gff3 into genomic coordinates
  transdecoder=coordfix(transdecoder,cds,transeq)
  transdecoder,trans=converter(transcripts,transdecoder,genome)
# TRANSDECODER ONLY
# Step 3. Print gtf, run seqret, load proteins, modify genbank
  gtfout(transdecoder,TRANSDECODEROUT+".gtf")
  `cat #{TRANSDECODEROUT}.gtf #{TRANSCRIPTS} | sort -k1,1 -k 4,4n |sed 's/gene_id/ID/'| ruby -ne 'line=$_.chomp.split; line.map!{|x|x.chomp(";").gsub(/"/,"")};puts(line[0..7].join("\t")+"\tlocus_tag="+line[8].split("=")[1].split(";")[0]+"; "+line[8..-1].join("; "))' | sed 's/transcript/gene/' > merged.gtf`
  `grep 'NC_003030.1' merged.gtf > #{FINALCHROM}`
  `grep 'NC_001988.2' merged.gtf > #{FINALPLAS}`
  seqret(CHROMFASTA,FINALCHROM)
  seqret(PLASFASTA,FINALPLAS)

  `cat nc_003030.genbank nc_001988.genbank > #{TEMP}`
  proteins=proteinextract(transdecoder,genome)
  readnwrite(TEMP,TRANSDECODEROUT+".genbank",transdecoder,proteins)
# Standard only (for revisions)
  `cat #{TRANSCRIPTS} #{STANDARDCDS} | sort -k1,1 -k 4,4n | sed 's/gene_id/ID/' | ruby -ne 'line=$_.chomp.split; line.map!{|x|x.chomp(";").gsub(/"/,"")}; out=line[0..7].join("\t")+"\tlocus_tag="; out+=(line[8].include?("=") ? line[8].split("=")[1]+"; "+line[8..-1].join("; ") : line[9]+"; ID="+line[9]); out+=";"; puts(out)' | sed 's/transcript/gene/' > merged.gtf`
  `grep 'NC_003030.1' merged.gtf > #{FINALCHROM}`
  `grep 'NC_001988.2' merged.gtf > #{FINALPLAS}`
  seqret(CHROMFASTA,FINALCHROM)
  seqret(PLASFASTA,FINALPLAS)
  `cp #{STANDARDCDS} #{STANDARDOUT}.gtf`
  `cat nc_003030.genbank nc_001988.genbank > #{TEMP}`
  #proteins=proteinextract(transdecoder,genome)
  readnwrite(TEMP,STANDARDOUT+".genbank",transdecoder,proteins)
# TOTAL
  `cat #{TRANSDECODEROUT}.gtf #{TRANSCRIPTS} #{STANDARDCDS} | sort -k1,1 -k 4,4n | sed 's/gene_id/ID/' | ruby -ne 'line=$_.chomp.split; line.map!{|x|x.chomp(";").gsub(/"/,"")}; out=line[0..7].join("\t")+"\tlocus_tag="; out+=(line[8].include?("=") ? line[8].split("=")[1]+"; "+line[8..-1].join("; ") : line[9]+"; ID="+line[9]); out+=";"; puts(out)' | sed 's/transcript/gene/' > merged.gtf`
  `grep 'NC_003030.1' merged.gtf > #{FINALCHROM}`
  `grep 'NC_001988.2' merged.gtf > #{FINALPLAS}`
  seqret(CHROMFASTA,FINALCHROM)
  seqret(PLASFASTA,FINALPLAS)
  `cat #{TRANSDECODEROUT}.gtf #{STANDARDCDS} > #{MERGEDOUT}.gtf`
  `cat nc_003030.genbank nc_001988.genbank > #{TEMP}`
  proteins=proteinextract(transdecoder,genome)
  readnwrite(TEMP,MERGEDOUT+".genbank",transdecoder,proteins)


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
