#!/bin/bash
#PBS -N visualization
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=15
#PBS -l walltime=80:00:00
#PBS -d /home/mrals/Final
#------------------------------------------------
# Title: viewing.sh
#
# Copyright 2014 Matt Ralston
# 
# This script assists in the visualization of
# coverage and gene expression. 
#
#                P A R T    1
# Here alignment files are split by strand,
# where paired and unpaired reads are grouped by
# ScriptSeq v2 strand-specificity. These files
# can be visualized in IGV, for example.
#
#                P A R T    2
# Here, alignment files are parsed into BED format
# files for read counting and visualization in
# circos.
# 
#------------------------------------------------


#set -e

#------------------------------------------------
# Source
#------------------------------------------------
source functions.sh
#------------------------------------------------
# Parameters
#------------------------------------------------

# INDIR: This is the location of the processed
# alignment files to be visualized.
INDIR=SAM_processed
# OUTDIR: This is where the strand-specific
# alignment files will be produced
export OUTDIR=BAM_strand
# MRG: This is the location where temporary
# alignment files will be created for merging.
export MRG="tmp/merge"
# REFGENOME: This is the location of a .genome file
# used for bedtools.
export REFGENOME=reference/CAC.genome
# CORES: This is the number of processor cores to
# be used for parallelization.
CORES=15
# EXPRDIR: This is where the bed and bedgraph coverage
# files will be produced
export EXPRDIR=Expression
# CIRC: This is the location where all of the circos
# plots and configuration files will be produced.
CIRC=/home/mrals/Final/circos
# BASE: This is the home directory of the script.
base=/home/mrals/Final
export base
# PICARD: This is the location of picard jar files.
export PICARD=/usr/local/picard-tools-1.67
# ANNOTATION: This is the location of a bed format
# CDS annotation. May be deprecated.
ANNOTATION=reference/CAC.bed
FILES=`/usr/bin/ls $INDIR/*.3.bam`

#################    P A R T   1   ###################

parallel -j $CORES 'bamsplit {} $MRG $PICARD $REFGENOME $OUTDIR $EXPRDIR' ::: $FILES


#./coverage_calc.rb



#samtools mpileup -BQ0 run.sorted.bam | perl -pe '($c, $start, undef, $depth) = split;if ($c ne $lastC || $start != $lastStart+1) {print "fixedStep chrom=$c start=$start step=1 span=1\n";}$_ = $depth."\n";($lastC, $lastStart) = ($c, $start);' | gzip -c > run.wig.gz





## EOF-------------------------------------------
