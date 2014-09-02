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


function bamsplit {
# This function takes a 'full' BAM file as a primary argument
# Merge directory, location of Picard .jars, reference genome,
# stranded output directory, and expression directory as secondary arguments
    file=$1
    MRG=$2
    PICARD=$3
    REFGENOME=$4
    OUTDIR=$5
    EXPRDIR=$6
    f=${file##*/};f=${f%*.*.*}
############################### P A I R E D ########################
    ########      P L U S
    # Discordant alignment type 1: dis1
    # this first type of alignment refers to a discordant type alignment
    # to the positive strand, where the first read aligns to the + strand
    # and the second read does not align to either strand (discordant)
    samtools view -hb -f 73 $file | samtools view -hb -F 20 - > $MRG/pair.plus.for.dis1
    samtools view -hb -f 133 $file | samtools view -hb -F 40 - > $MRG/pair.plus.rev.dis1
    # Discordant alignment type 2: dis 2
    # The second type of discordant alignment is where the first read does not
    # align to either strand and the second read aligns to the minus strand.
    samtools view -hb -f 101 $file | samtools view -hb -F 8 - > $MRG/pair.plus.for.dis2
    samtools view -hb -f 153 $file | samtools view -hb -F 36 - > $MRG/pair.plus.rev.dis2
    # Concordant alignment
    # This type of alignment is where the first read aligns to the forward strand
    # and the second read aligns to ther second strand, indicating that the sequenced
    # fragment originates from the forward (+) strand.
    samtools view -hb -f 97 $file | samtools view -hb -F 28 - > $MRG/pair.plus.for.con
    samtools view -hb -f 145 $file | samtools view -hb -F 44 - > $MRG/pair.plus.rev.con
    
    ########      M I N U S
    # Discordant alignment type 1: dis1
    # this first type of alignment refers to a discordant type alignment
    # to the negative strand, where the first read aligns to the - strand
    # and the second read does not align to either strand (discordant)

    samtools view -hb -f 89 $file | samtools view -hb -F 36 - > $MRG/pair.minus.for.dis1
    samtools view -hb -f 165 $file | samtools view -hb -F 24 - > $MRG/pair.minus.rev.dis1
    # Discordant alignment type 2: dis 2
    # The second type of discordant alignment is where the first read does not
    # align to either strand and the second read aligns to the + strand.
    samtools view -hb -f 69 $file | samtools view -hb -F 40 - > $MRG/pair.minus.for.dis2
    samtools view -hb -f 137 $file | samtools view -hb -F 52 - > $MRG/pair.minus.rev.dis2
    # Concordant alignment
    # This type of alignment is where the first read aligns to the - strand
    # and the second read aligns to the + strand, indicating that the sequenced
    # fragment originates from the reverse (-) strand.
    samtools view -hb -f 81 $file | samtools view -hb -F 44 - > $MRG/pair.minus.for.con
    samtools view -hb -f 161 $file | samtools view -hb -F 28 - > $MRG/pair.minus.rev.con
    ############# U N P A I R E D #########
    samtools view -hb -F 21 $file > $MRG/un.plus
    samtools view -hb -F 5 $file | samtools view -hb -f 16 - > $MRG/un.minus
    # Merge - sort
    plus=`/usr/bin/ls $MRG/*.plus*`
    minus=`/usr/bin/ls $MRG/*.minus*`
    samtools merge -f tmp/plus.bam $plus 
    samtools sort tmp/plus.bam $OUTDIR/$f.plus > $OUTDIR/$f.plus
    samtools merge -f tmp/minus.bam $minus
    samtools sort tmp/minus.bam $OUTDIR/$f.minus > $OUTDIR/$f.minus
    rm tmp/plus.bam tmp/minus.bam
    samtools index $OUTDIR/$f.plus.bam $OUTDIR/$f.plus.bam.bai TMP_DIR=tmp
    samtools index $OUTDIR/$f.minus.bam $OUTDIR/$f.minus.bam.bai TMP_DIR=tmp
    #      P R O C E S S   B A M   F O R  C O V E R A G E
    #          B E D G R A P H
    java -jar $PICARD/MarkDuplicates.jar INPUT=$OUTDIR/$f.plus.bam OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | samtools view -b -h -q 20 - | bedtools genomecov -bg -ibam /dev/stdin -g $REFGENOME > $EXPRDIR/$f+.bedgraph
    java -jar $PICARD/MarkDuplicates.jar INPUT=$OUTDIR/$f.minus.bam OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | samtools view -b -h -q 20 - | bedtools genomecov -bg -ibam /dev/stdin -g $REFGENOME > $EXPRDIR/$f-.bedgraph
    #          C O V E R A G E
    # option to set genomecov option to -bg (areas of continuous coverage together)
    java -jar $PICARD/MarkDuplicates.jar INPUT=$OUTDIR/$f.plus.bam OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | samtools view -b -h -q 20 - | bedtools genomecov -d -ibam /dev/stdin -g $REFGENOME > $EXPRDIR/$f+.cov
    java -jar $PICARD/MarkDuplicates.jar INPUT=$OUTDIR/$f.minus.bam OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | samtools view -b -h -q 20 - | bedtools genomecov -d -ibam /dev/stdin -g $REFGENOME > $EXPRDIR/$f-.cov
    sed -i "s/NC_003030.1/C/g" $EXPRDIR/$f+.cov
    sed -i "s/NC_001988.2/P/g" $EXPRDIR/$f+.cov
    sed -i "s/NC_003030.1/C/g" $EXPRDIR/$f-.cov
    sed -i "s/NC_001988.2/P/g" $EXPRDIR/$f-.cov
    cp $EXPRDIR/$f+.cov $EXPRDIR/$f+.log.cov
    cp $EXPRDIR/$f-.cov $EXPRDIR/$f-.log.cov
    ruby -p -i -e 'temp=$_.chomp.split; puts (temp[0..2].join("\t")+"\t#{Math.log10(temp[3].to_i+1)}")' $EXPRDIR/$f+log.bedgraph
    ruby -p -i -e 'temp=$_.chomp.split; puts (temp[0..2].join("\t")+"\t#{Math.log10(temp[3].to_i+1)}")' $EXPRDIR/$f-log.bedgraph

}



export -f quality preprocess bamsplit

## EOF------------------------------------------
