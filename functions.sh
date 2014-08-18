#!/bin/bash

function quality {
# This function processes the quality of fastq.gz files with FastQC, Prinseq, and the fastx toolkit

# This function takes a file as a primary argument,
# an input directory as a secondary argument,
# and a output directory as a tertiary argument
    file=$1
    indir=$2
    outdir=$3
    mkdir $outdir/${file%.*}
    fastqc -j /usr/bin/java -f fastq -o $outdir/${file%.*} $indir/${file%.*}
    zcat $indir/$file | fastx_quality_stats -Q 33 > $outdir/${file%.*}/fastx_report.txt
    zcat $indir/$file | prinseq -fastq stdin -stats_all > $outdir/${file%.*}/prinseq_stats.txt
}

function preprocess {
# This function preprocesses fastq.gz files to concatenate the Casava 1.8+ header
# Which is normally split by whitespace:
# @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG

# This function takes a file as a primary argument
# and a input directory as a secondary argument.
    file=$1
    indir=$2
    gunzip $indir/$file
    cat $indir/${file%.*} | ruby -ne 'puts($_.chomp.split.join("_"))' | gzip > $indir/$file
    rm $indir/${file%.*}
}



export -f quality preprocess

## EOF------------------------------------------
