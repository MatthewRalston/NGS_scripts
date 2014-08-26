#!/bin/bash
#PBS -N trinity
#PBS -r n
#PBS -V
#PBS -l nodes=biohen29:ppn=20
#PBS -l walltime=100:00:00
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

#------------------------------------------------
# Parameters
#------------------------------------------------
# CORES
CORES=20
JM=30G
# Trinity

# Reference genome
REFGENOME=reference/CAC.txt
#
BAMDIR=SAM_processed
INDIR=/home/mrals/Final/final/finalfastq
OUTDIR=/home/mrals/Final/Trinity_ref
#OUTDIR=/home/mrals/Final/Trinity_de_novo
#TRIN=Trinity.fasta
TRIN=Trinity-GG.fasta
TMP=tmp
RAWDIR=rawdata
FILES1=`/usr/bin/ls $RAWDIR/*_1.fastq`
FILES2=`/usr/bin/ls $RAWDIR/*_2.fastq`
R1='$_[0..3] == "@HWI" ? puts("#{$_.chomp/1") : puts($_.chomp)'
R2='$_[0..3] == "@HWI" ? puts("#{$_.chomp/2") : puts($_.chomp)'
export R1 R2

##################################################
#
# Merge fastq files
#
##################################################
#zcat $INDIR/*.unpaired.gz | grep " 1:N:0" | gzip >> $OUTDIR/t1.fq
#zcat $INDIR/*.unpaired.gz | grep " 2:N:0" | gzip >> $OUTDIR/t2.fq
#zcat $OUTDIR/t1.fq $INDIR/*.1.gz | ruby -ne '$_[0..2] == "@HWI" ? puts("#{$_.chomp}/1") : puts($_)' | gzip > $OUTDIR/left.fq.gz
#zcat $OUTDIR/t2.fq $INDIR/*.2.gz | ruby -ne '$_[0..2] == "@HWI" ? puts("#{$_.chomp}/2") : puts($_)' | gzip > $OUTDIR/right.fq.gz
#rm $OUTDIR/t1.fq $OUTDIR/t2.fq


##################################################
#
# OPTIONAL-- Reference Trinity
# Use BAM files (instead of fastq files
##################################################
samtools merge All.merged.bam `/usr/bin/ls -d $BAMDIR/*.3.bam` 
samtools sort -o $TMP/All.merged.bam $TMP/All.most
# stupid b.s. to add the g.d. suffixes.
samtools view -h $TMP/All.most | ruby -ne '$_[0..2] == "HWI" ? puts($_.split[0]+"/"+$_.split("_")[1][0]+"\t"+$_.split[1..$_.split.size].join("\t")) : puts($_.chomp)' | samtools view -bh - $TMP/All.bam
rm All.merged.bam


###Trinity --genome $REFGENOME --genome_guided_use_bam $TMP/All.bam --genome_guided_max_intron 1 --genome_guided_sort_buffer 15G --genome_guided_CPU $CORES --SS_lib_type FR --seqType fq --jaccard_clip --JM $JM --CPU $CORES --output $OUTDIR --left $TMP/left.fq --right $TMP/right.fq &> $OUTDIR/ref_assembly.log

#Trinity --left $TMP/left.fq --right $TMP/right.fq --jaccard_clip --genome $REFGENOME --genome_guided_max_intron 1 --genome_guided_use_bam $TMP/FINAL.bam --JM $JM --seqType fq --output $OUTDIR --genome_guided_CPU $CORES --CPU $CORES --SS_lib_type FR


##################################################
#
# De-novo Trinity
# This step performs the initial assembly with trinity, producing a fasta assembly Trinity.fa
# in $OUTDIR
##################################################
#Trinity --seqType fq --JM $JM --left $OUTDIR/left.fq --right $OUTDIR/right.fq --SS_lib_type FR --CPU $CORES --jaccard_clip --output $OUTDIR &> $OUTDIR/denovo_assembly.log

##################################################
#
# FASTA -> GTF
# This step transforms this assembly into a gtf format assembly, although the features/CDSes
# of the traditional CAC genes are not mapped back on to these features yet...
##################################################
#bwa mem -t $CORES $REFGENOME $OUTDIR/$TRIN | samtools view -Sbh - | samtools sort - $OUTDIR/Trinity
# Indexing the alignment
#samtools index $OUTDIR/Trinity.bam
# Conversion to gtf, gene_ids will be Trinity transcripts, though
#samtools view -bh $OUTDIR/Trinity.bam |  bam2bed | ruby -ne 'puts($_.split("\t")[0...6].join("\t"))' | bedToGenePred stdin stdout| genePredToGtf file stdin $OUTDIR/Trinity.gtf

##################################################
#
# Annotation merging
# This step uses a custom ruby script to merge the assembly and the CDS annotation
##################################################






##################################################
#
# Summary
# 
##################################################
#transrate



## EOF-------------------------------------------
