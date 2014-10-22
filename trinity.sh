#!/bin/env bash
#PBS -N trinity
#PBS -r n
#PBS -V
#PBS -l nodes=biohen29:ppn=20
#PBS -l walltime=100:00:00
#PBS -l mem=120gb
#PBS -d /home/mrals/Final
#------------------------------------------------
# Title: trinity.sh
#
# Copyright 2014 Matt Ralston
#
# This script performs a de-novo transcriptome
# assembly and converts this into gtf format.
# It turns out that trinity can be rather difficult to use.
# The simplest way to get it to work is to go back to the
# Illumina Fastq file, and concatenate the read number and index
# i.e. the 1:N:0:TTAGGC type suffix to the read location with 
# an underscore: _  Also, bowtie2 will remove any trailing
# indicators for the reads /1, so it may be necessary to append this.
#
#------------------------------------------------


#set -e

source functions.sh
#------------------------------------------------
# Parameters
#------------------------------------------------
# CORES
export CORES=20
export JM=5G
# Trinity
WORKDIR=/home/mrals/Final
# Reference genome
export REFGENOME=$WORKDIR/reference/CAC.txt
# Reference proteome
export REFPROTEOME=$WORKDIR/reference/CAC_proteins.fasta
#
export BAMDIR=$WORKDIR/SAM_processed
export INDIR=$WORKDIR/final/finalfastq
#REFOUT=$WORKDIR/Trinity_ref
#REFOUT=$WORKDIR/Trinity_paired
#REFOUT=$WORKDIR/Trinity_de_novo
#TRIN=Trinity.fasta
RAW=Trin_raw
export TRIN=Trinity-GG.fasta
TRINDN=Trinity.fasta
TMP=tmp
#RAWDIR=rawdata
#FILES1=`/usr/bin/ls $RAWDIR/*_1.fastq.gz`
#FILES2=`/usr/bin/ls $RAWDIR/*_2.fastq.gz`
R1='$_[0..3] == "@HWI" ? puts("#{$_.chomp/1") : puts($_.chomp)'
R2='$_[0..3] == "@HWI" ? puts("#{$_.chomp/2") : puts($_.chomp)'
COMPARE=("NS*15" "NS*75" "NS*150" "NS*270" "BuOH*15" "BuOH*75" "BuOH*150" "BuOH*270" "BA*15" "BA*75" "BA*150" "BA*270")

export R1 R2

##################################################
#
# Individual
#
##################################################



#for files in ${COMPARE[@]};
#do
#    cd Assemblies
#    mintrin $files $INDIR $BAMDIR $REFGENOME
#    cd ..
#done

##################################################
#
# Merge fastq files
#
##################################################
#zcat $INDIR/*unpaired.gz | ./untangler.rb 1 | gzip > $RAW/t1.fq.gz
#zcat $INDIR/*unpaired.gz | ./untangler.rb 2 | gzip > $RAW/t2.fq.gz
# Transrate only
#zcat $INDIR/*unpaired.gz | ruby -ne '$_[0..3] == "@HWI" ? puts($_.chomp.split("_")[0]) : puts($_.chomp)' | gzip > $RAW/unpaired.fq.gz
#zcat $INDIR/*.1.gz | ruby -ne '$_[0..3] == "@HWI" ? puts($_.chomp.split("_")[0]) : puts($_.chomp)' | gzip > $RAW/left.fq.gz
#zcat $INDIR/*.2.gz | ruby -ne '$_[0..3] == "@HWI" ? puts($_.chomp.split("_")[0]) : puts($_.chomp)' | gzip > $RAW/right.fq.gz
# Trinity only
#zcat $INDIR/*.1.gz $RAW/t1.fq.gz | ruby -ne '$_[0..3] == "@HWI" ? puts($_.chomp+"/1") : puts($_.chomp)' | gzip > $RAW/left_combined.fq.gz
#zcat $INDIR/*.2.gz | ruby -ne '$_[0..3] == "@HWI" ? puts($_.chomp+"/2") : puts($_.chomp)' | gzip > $RAW/right_combined.fq.gz
# Corrects an issue with all unpaired right reads being reverse complemented in trimfiltercheck
#zcat $RAW/t2.fq.gz | ruby -ne '$_[0..3] == "@HWI" ? puts($_.chomp+"/2") : puts($_.chomp)' | fastx_reverse_complement -Q33 | gzip >> $RAW/right_combined.fq.gz
#rm $RAW/t1.fq.gz $RAW/t2.fq.gz
##################################################
#
# OPTIONAL-- Reference Trinity
# Use BAM files (instead of fastq files
##################################################
# ALL DATA
#samtools merge $TMP/All.merged.bam `/usr/bin/ls -d $BAMDIR/*.3.bam` 
#samtools sort -o $TMP/All.merged.bam $TMP/All.most > $TMP/All.most.bam
#rm $TMP/All.merged.bam
# stupid b.s. to add the g.d. suffixes.
#samtools view -h $TMP/All.most.bam | ruby -ne '$_[0..2] == "HWI" ? puts($_.split[0]+"/"+$_.split("_")[1][0]+"\t"+$_.split[1..$_.split.size].join("\t")) : puts($_.chomp)' | ruby -ne '$_[0..2] == "HWI" ? (  l=$_.split("\t"); l[1].to_i%2 == 1 ? puts(l.join("\t")) : puts( ([l[0]]+[(l[1].to_i+1).to_s]+l[2..l.size]).join("\t") )  ) : puts($_.chomp)' | samtools view -hbS - > $RAW/All.bam
#Trinity --left $RAW/left_combined.fq.gz --right $RAW/right_combined.fq.gz --genome $REFGENOME --genome_guided_max_intron 1 --jaccard_clip --genome_guided_use_bam $RAW/All.bam --JM $JM --seqType fq --output $REFOUT --genome_guided_CPU 4 --CPU $CORES --SS_lib_type FR >> $REFOUT/ref_assembly.log 2>&1

# PAIRED ONLY
#samtools view -h $TMP/All.most.bam | ruby -ne '$_[0..2] == "HWI" ? puts($_.split[0]+"/"+$_.split("_")[1][0]+"\t"+$_.split[1..$_.split.size].join("\t")) : puts($_.chomp)' | ruby -ne '$_[0..2] == "HWI" ? (  l=$_.split("\t"); l[1].to_i%2 == 1 ? puts(l.join("\t")) : puts( ([l[0]]+[(l[1].to_i+1).to_s]+l[2..l.size]).join("\t") )  ) : puts($_.chomp)' | samtools view -hbS -f 2 - > $RAW/All.bam
#Trinity --left $RAW/left.fq.gz --right $RAW/right.fq.gz --genome $REFGENOME --genome_guided_max_intron 1 --jaccard_clip --genome_guided_use_bam $RAW/All.bam --JM $JM --seqType fq --output $REFOUT --genome_guided_CPU 4 --CPU $CORES --SS_lib_type FR >> $REFOUT/ref_assembly.log 2>&1


##################################################
#
# De-novo Trinity
# This step performs the initial assembly with trinity, producing a fasta assembly Trinity.fa
# in $OUTDIR
##################################################
#Trinity --seqType fq --JM $JM --left $RAW/left_combined.fq.gz --right $RAW/right_combined.fq.gz --SS_lib_type FR --CPU $CORES --jaccard_clip --output $TRINOUT &> $TRINOUT/denovo_assembly.log
##################################################
#
#    A s s e m b l y    M e t r i c s
#
##################################################
#  R e f e r e n c e
#rvm use ruby-2.1.2@transrate
#transrate -a transrate/$TRIN -r $REFPROTEOME -g $REFGENOME -l $RAW/left.fq.gz -i $RAW/right.fq.gz -u $RAW/unpaired.fq.gz -s fr -o transrate/singletons.sam -f transrate/transrate_output.csv -t $CORES -x 0
# QUICK
#transrate -a transrate/$TRIN -r $REFPROTEOME -g $REFGENOME -x 1 -f transrate/transrate_output.csv
#  D e    N o v o
#transrate -a $TRINOUT/$TRINDN -r $REFPROTEOME -g $REFGENOME -x 1 -f $TRINOUT/transrate_output.csv
##################################################
#
# FASTA -> GTF
# This step transforms this assembly into a gtf format assembly, although the features/CDSes
# of the traditional CAC genes are not mapped back on to these features yet...
##################################################
#sed -i 's/|/-/' $REFOUT/$TRIN


#bwa mem -t $CORES $REFGENOME $REFOUT/$TRIN | samtools view -Sbh - | samtools sort - $REFOUT/Trinity
# Indexing the alignment
#samtools index $REFOUT/Trinity.bam
# Conversion to gtf, gene_ids will be Trinity transcripts, though
#samtools view -bh $REFOUT/Trinity.bam |  bam2bed | ruby -ne 'puts($_.split("\t")[0...6].join("\t"))' | bedToGenePred stdin stdout| genePredToGtf file stdin $REFOUT/Trin-bwa.gtf


blat $REFGENOME $REFOUT/$TRIN -maxIntron=0 $REFOUT/Trinity.psl
psl2bed-best $REFOUT/Trinity.psl /dev/stdout | ruby -ne 'puts($_.split("\t")[0...6].join("\t"))' | bedToGenePred stdin stdout| genePredToGtf file stdin $REFOUT/Trin-blat.gtf


blastn -query $REFOUT/$TRIN -db $REFGENOME -reward 2 -penalty -3 -gapopen 6 -gapextend 2 -outfmt 6 > $REFOUT/Trinity-GG.blast
./sortblast.rb
# The curate partial and combine with unique transcripts
#cat $REFOUT/transcripts-unique.gtf | ./assemblycurate.rb | sort -k 1,1 -k4,4n > $REFOUT/Trin-blast.gtf



#blat $REFGENOME $TRINOUT/$TRINDN -maxIntron=0 $TRINOUT/Trinity.psl
#psl2bed-best $TRINOUT/Trinity.psl /dev/stdout | ruby -ne 'puts($_.split("\t")[0...6].join("\t"))' | bedToGenePred stdin stdout| genePredToGtf file stdin $TRINOUT/Trinity.gtf
##################################################
#
# GTF merging
# This step uses a custom ruby script to merge the assembly and the CDS annotation
##################################################
#cat reference/CAC.gtf | grep 'CA_Ct' > reference/CAC.trna.gtf
#cat reference/CAC.gtf | grep 'CA_Cr' > reference/CAC.rrna.gtf
#cat reference/CAC.gtf | grep -v 'CA_Cr' | grep -v 'CA_Ct' | grep 'exon' | sed 's/exon/CDS/' > reference/CAC.CDS.gtf
#      R E F E R E N C E
#cat $REFOUT/Trinity.gtf | grep -v 'codon' | grep -v 'CDS' | sed 's/stdin/trinity/' | sed 's/exon\t/transcript\t/' | sed 's/|/-/' | ./assemblycurate.rb | sort -k 1,1 -k4,4n > $REFOUT/Trinity_filtered.gtf
#cat $REFOUT/Trinity_filtered.gtf reference/CAC.trna.gtf reference/CAC.rrna.gtf reference/CAC.CDS.gtf | sort -k 1,1 -k 4,4n > $REFOUT/combined.gtf
#      D E    N O V O
#cat $TRINOUT/Trinity.gtf | grep -v 'codon' | grep -v 'CDS' | sed 's/stdin/trinity/' | sed 's/exon\t/transcript\t/' | ./assemblycurate.rb > $TRINOUT/Trinity_filtered.gtf
#cat $REFOUT/Trinity_filtered.gtf reference/CAC.trna.gtf reference/CAC.rrna.gtf reference/CAC.CDS.gtf | sort -k 1,1 -k 4,4n > $TRINOUT/combined.gtf





## EOF-------------------------------------------
