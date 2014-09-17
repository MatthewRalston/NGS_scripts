#!/bin/bash
#PBS -N trinity
#PBS -r n
#PBS -V
#PBS -l nodes=biohen29:ppn=20
#PBS -l walltime=100:00:00
#PBS -l mem=100gb
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
JM=4G
# Trinity

# Reference genome
REFGENOME=reference/CAC.txt
# Reference proteome
REFPROTEOME=reference/CAC_proteins.fasta
#
BAMDIR=SAM_processed
INDIR=/home/mrals/Final/final/finalfastq
OUTDIR=/home/mrals/Final/Trinity_ref
#OUTDIR=/home/mrals/Final/Trinity_de_novo
#TRIN=Trinity.fasta
RAW=Trin_raw
TRIN=Trinity-GG.fasta
TMP=tmp
#RAWDIR=rawdata
#FILES1=`/usr/bin/ls $RAWDIR/*_1.fastq.gz`
#FILES2=`/usr/bin/ls $RAWDIR/*_2.fastq.gz`
R1='$_[0..3] == "@HWI" ? puts("#{$_.chomp/1") : puts($_.chomp)'
R2='$_[0..3] == "@HWI" ? puts("#{$_.chomp/2") : puts($_.chomp)'
export R1 R2

##################################################
#
# Merge fastq files
#
##################################################
#zcat $INDIR/*unpaired.gz | grep "1:N:0" | gzip >> $OUTDIR/t1.fq.gz
#zcat $INDIR/*unpaired.gz | grep "2:N:0" | gzip >> $OUTDIR/t2.fq.gz
#cat $OUTDIR/t*.fq.gz  > $RAW/unpaired.fq.gz
#cat $INDIR/*.1.gz > $RAW/left.fq.gz
#cat $INDIR/*.2.gz > $RAW/right.fq.gz
#zcat $OUTDIR/t1.fq.gz $INDIR/*.1.gz | ruby -ne '$_[0..2] == "@HWI" ? puts("#{$_.chomp}/1") : puts($_)' | gzip > $RAW/left_combined.fq.gz
#zcat $OUTDIR/t2.fq.gz $INDIR/*.2.gz | ruby -ne '$_[0..2] == "@HWI" ? puts("#{$_.chomp}/2") : puts($_)' | gzip > $RAW/right_combined.fq.gz
#rm $OUTDIR/t1.fq.gz $OUTDIR/t2.fq.gz


##################################################
#
# OPTIONAL-- Reference Trinity
# Use BAM files (instead of fastq files
##################################################
#samtools merge $TMP/All.merged.bam `/usr/bin/ls -d $BAMDIR/*.3.bam` 
#samtools sort -o $TMP/All.merged.bam $TMP/All.most > $TMP/All.most.bam
#rm $TMP/All.merged.bam
# stupid b.s. to add the g.d. suffixes.
#samtools view -h $TMP/All.most.bam | ruby -ne '$_[0..2] == "HWI" ? puts($_.split[0]+"/"+$_.split("_")[1][0]+"\t"+$_.split[1..$_.split.size].join("\t")) : puts($_.chomp)' | ruby -ne '$_[0..2] == "HWI" ? (  l=$_.split("\t"); l[1].to_i%2 == 1 ? puts(l.join("\t")) : puts( ([l[0]]+[(l[1].to_i+1).to_s]+l[2..l.size]).join("\t") )  ) : puts($_.chomp)' | samtools view -hbS - > $OUTDIR/All.bam


###Trinity --genome $REFGENOME --genome_guided_use_bam $TMP/All.bam --genome_guided_max_intron 1 --genome_guided_sort_buffer 15G --genome_guided_CPU $CORES --SS_lib_type FR --seqType fq --jaccard_clip --JM $JM --CPU $CORES --output $OUTDIR --left $TMP/left.fq --right $TMP/right.fq &> $OUTDIR/ref_assembly.log

#Trinity --left $RAW/left_combined.fq.gz --right $RAW/right_combined.fq.gz --jaccard_clip --genome $REFGENOME --genome_guided_max_intron 1 --genome_guided_use_bam $RAW/All.bam --JM $JM --seqType fq --output $OUTDIR --genome_guided_CPU 4 --CPU $CORES --SS_lib_type FR


##################################################
#
# De-novo Trinity
# This step performs the initial assembly with trinity, producing a fasta assembly Trinity.fa
# in $OUTDIR
##################################################
#Trinity --seqType fq --JM $JM --left $OUTDIR/left.fq --right $OUTDIR/right.fq --SS_lib_type FR --CPU $CORES --jaccard_clip --output $OUTDIR &> $OUTDIR/denovo_assembly.log


##################################################
#
#    A s s e m b l y    M e t r i c s
#
##################################################
transrate -a $OUTDIR/$TRIN -r $REFPROTEOME -g $REFGENOME -l $RAW/left.fq.gz -i $RAW/right.fq.gz -u $RAW/unpaired.fq.gz -s fr -o $OUTDIR/singletons.sam -f $OUTDIR/transrate_output.csv -t $CORES -x 0







##################################################
#
# FASTA -> GTF
# This step transforms this assembly into a gtf format assembly, although the features/CDSes
# of the traditional CAC genes are not mapped back on to these features yet...
##################################################

#blat $REFGENOME $OUTDIR/$TRIN -maxIntron=0 $OUTDIR/Trinity.psl
#psl2bed --headered < $OUTDIR/Trinity.psl | ruby -ne 'puts($_.split("\t")[0...6].join("\t"))' | bedToGenePred stdin stdout| genePredToGtf file stdin $OUTDIR/Trinity.gtf

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




## EOF-------------------------------------------
