#!/bin/bash
#PBS -N trinity
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=6
#PBS -l walltime=240:00:00
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
CORES=6
JM=30G
# Trinity

# Reference genome
REFGENOME=reference/CAC.txt
#
BAMDIR=SAM_processed
INDIR=/home/mrals/Final/final/finalfastq
OUTDIR=/home/mrals/Final/Trinity
#OUTDIR=/home/mrals/Final/ggtrin
TRIN=Trinity.fasta
#TRIN=Trinity-GG.fasta
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
cat $INDIR/*.1.gz > $OUTDIR/t1.fq
cat $INDIR/*.2.gz > $OUTDIR/t2.fq
cat $INDIR/*.unpaired.gz | grep " 1:N:0" >> $OUTDIR/t1.fq
cat $INDIR/*.unpaired.gz | grep " 2:N:0" >> $OUTDIR/t2.fq
cat $OUTDIR/t1.fq | ruby -ne '$_[0..2] == "@HWI" ? puts("#{$_.chomp}/1") : puts($_)' >  $OUTDIR/left.fq
cat $OUTDIR/t2.fq | ruby -ne '$_[0..2] == "@HWI" ? puts("#{$_.chomp}/2") : puts($_)' >  $OUTDIR/right.fq
rm $OUTDIR/t1.fq $OUTDIR/t2.fq
# NOTE Fastq files must have /1 or /2 preappendend to identifiers where appropriate
#parallel -j $CORES "cat {} | ruby -ne '$R1' | gzip > {}.gz" ::: $FILES1
#parallel -j $CORES "cat {} | ruby -ne '$R2' | gzip > {}.gz" ::: $FILES2


##################################################
#
# OPTIONAL-- Reference Trinity
# Use BAM files (instead of fastq files
##################################################
samtools merge All.merged.bam `/usr/bin/ls -d $BAMDIR/*.3.bam` 
samtools sort -o $TMP/All.merged.bam $TMP/All.most
# stupid b.s. to add the g.d. suffixes.
samtools view -h $TMP/All.bam | ruby -ne '$_[0..2] == "HWI" ? puts($_.split[0]+"/"+$_.split("_")[1][0]+"\t"+$_.split[1..$_.split.size].join("\t")) : puts($_.chomp)' | samtools view -bh - $TMP/All.bam


#Trinity --genome $REFGENOME --genome_guided_use_bam $TMP/All.bam --genome_guided_max_intron 1 --genome_guided_sort_buffer 15G --genome_guided_CPU $CORES --SS_lib_type FR --seqType fq --jaccard_clip --JM $JM --CPU $CORES --output $OUTDIR --left $TMP/left.fq --right $TMP/right.fq &> $OUTDIR/ref_assembly.log
Trinity --left $TMP/left.fq --right $TMP/right.fq --jaccard_clip --genome $REFGENOME --genome_guided_max_intron 1 --genome_guided_use_bam $TMP/FINAL.bam --JM $JM --seqType fq --output $OUTDIR --genome_guided_CPU $CORES --CPU $CORES --SS_lib_type FR


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
# Singletons
# The resulting report ~/ETP/trinity.oXXXXX will contain the information about the singletons
##################################################
#bowtie2-build $OUTDIR/$TRIN Trinity
#bowtie2 -p $CORES --fast -x Trinity --un-gz  $OUTDIR/singletons -1 left.fq -2 right.fq -S /dev/null
#tmp=$(zcat $OUTDIR/singletons. | wc -l )
#echo $(($tmp / 4)) > $OUTDIR/singletons.txt

##################################################
#
# Summary
# 
##################################################




## EOF-------------------------------------------
