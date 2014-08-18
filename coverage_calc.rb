#!/home/mrals/.rvm/rubies/ruby-2.0.0-p247/bin/ruby
=begin
------------------------------------------------------------------------------
--                                                                          --
--                                 MATTHEW RALSTON                          --
--                                                                          --
--                      C O V E R A G E _ C A L C. R B                      --
--                                                                          --
------------------------------------------------------------------------------
-- TITLE                                                                    --
--                                                                          --
--  Copyright 2014 Matthew Ralston                                          --
--                                                                          --
------------------------------------------------------------------------------
-- This file is designed to calculate coverage per gene, expressed as a     --
-- percentage, provided a base-by-base coverage (e.g. from BEDtools) file   --
-- and a gtf annotation from Cufflinks.					    --
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


#ANNOTATION="summary/summary2000.gtf"
ANNOTATION="reference/CAC.gtf"
DIRECTORY="Expression"
FASTA="reference/CAC.txt"
OUTDIR="summary/coverage/"
WINDOW=100

################################################
#
#               S U B - R O U T I N E S
#
################################################

################################################

################################################


def annotate
  genes={}
  File.open(ANNOTATION,'r').each do |line|
    temp=line.split
    genes[temp[9]]={}
    genes[temp[9]]["start"]=temp[3].to_i
    genes[temp[9]]["end"]=temp[4].to_i
    genes[temp[9]]["strand"]=temp[6]
    genes[temp[9]]["length"]=temp[4].to_i - 1 - temp[3].to_i
  end
  return genes
end

def avg(arr)
  ("%.2f" % (arr.inject {|sum,element| sum+element}.to_f / arr.size)).to_f
end

def stdev(arr)
  ("%.2f" % (Math.sqrt(avg(arr.map{|x| (x-avg(arr))**2})))).to_f
end

def genefinder(c, fastas, dict)
  gene=fastas[c][dict["start"]-1..dict["end"]]
  gene=gene.complement if dict["strand"] == "-"
  gene
end

def gccalc(gene)
  gene=Hash[*gene.split("").group_by{|i| i}.map{|k,v| [k,v.count]}.flatten]
  gene["c"] + gene["g"]
end


def percent_printer(avgcov,stdcov)
  avgcov.each do |cond,arr|
    File.open(OUTDIR+cond+".avcov",'w') do |file|
      file.puts("Gene_id\tlength\tgc\t"+(1..100).to_a.join("\t"))
      arr.each do |liszt|
        file.puts(liszt.join("\t"))
      end
    end
  end
  stdcov.each do |cond,arr|
    File.open(OUTDIR+cond+".sdcov",'w') do |file|
      arr.each do |liszt|
        file.puts(liszt.join("\t"))
      end
    end
  end
end


def tss_printer(tss)
  tss.each do |cond,arr|
    File.open(OUTDIR+cond+".tss",'w') do |file|
      file.puts("Gene_id\t"+(1...100).to_a.join("\t"))
      arr.each {|liszt| file.puts(liszt.join("\t"))}
    end
  end
end


def cov_calc(gcov)
	# calculate the coverage over the percentiles of the gene, given an array of per-base coverage for the gene.
  p=(1...100).to_a
  avcov=[];sdcov=[]
  # iterates through the percentiles, gathering the average and standard deviation of coverage throughout the percentile
  
  # avcov is a 100 element array : [ avg1, avg2, avg 3, ... etc ]
  p.each do |x|
    min=(gcov.size*(x-1)/100).floor
    max=(gcov.size*x/100).floor
    avcov << avg(gcov[min...max]); sdcov << stdev(gcov[min...max])
  end
  return avcov,sdcov
end


def over_genes(coverage,genes,strand,chrom,plas,fastas)
  # iterates through each gene, testing if the strand is appropriate
  
  # av is a 2D array [ 
  # gene1                   [ gene1, avg1, avg2, avg 3, ... etc ]
  # gene2                   [ gene2, avg1, avg2, avg 3, ... etc ]
  #                                     .      .      . 
  # geneX                   [ geneX, avg1, avg2, avg 3, ... etc ]
  #                                                            ]
  av=[];sd=[]; tsses=[]
  genes.each do |gene,dict|
    next if dict["strand"] != strand
    gene.include?("CA_C") ? x=chrom : x=plas
    gc=("%.2f" % (gccalc(genefinder(x,fastas,dict)).to_f/dict["length"]))
    # locates the coverage array index of the start and end of the gene
    gstart=(0...coverage.size).bsearch {|c| c >= dict["start"]}
    gend=(0...coverage.size).bsearch {|c| c >= dict["end"]+1}
    # initializes a subset of the coverage array, corresponding to the gene
    t1,t2=cov_calc(coverage[gstart...gend])
    av << [gene[1...-2], dict["length"], gc] + t1; sd << [gene[1...-2]]+ t2
    tsses << [gene[1...-2]] + coverage[gstart...gstart + WINDOW]
  end
  return av,sd,tsses
end


def over_files(fastas,genes,chrom,plas)
  avgcov={};stdcov={};tss_cov={}
  Dir.glob(DIRECTORY+"/*.cov") do |file|
    # initializes the coverage array for this file
    coverage=[]
    # tests for the strand of the coverage file.    
    file.include?("+") ? strand="+" : strand="-"
    # iterates through the file, including the coverage in the proper position of the array
    File.open(file,'r').each do |line|
      coverage << line.chomp.split[-1].to_i
    end
    avgcov[file.split("/")[-1].split(".")[0]], stdcov[file.split("/")[-1].split(".")[0]], tss_cov[file.split("/")[-1].split(".")[0]] = over_genes(coverage,genes,strand,chrom,plas,fastas)
  end
  percent_printer(avgcov,stdcov)
  tss_printer(tss_cov)
end


def main
	genes=annotate
  fastas={};  Bio::FastaFormat.open(FASTA).each_entry {|f| fastas[f.entry_id]=Bio::Sequence::NA.new(f.seq)}; chrom,plas = fastas.keys
  over_files(fastas,genes,chrom,plas)
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
