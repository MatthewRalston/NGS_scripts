#!/bin/bash
#PBS -N visualization
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=6
#PBS -l walltime=240:00:00
#PBS -d /home/mrals/ETP
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
# Parameters
#------------------------------------------------

# INDIR: This is the location of the processed
# alignment files to be visualized.
INDIR=SAM_processed
# OUTDIR: This is where the strand-specific
# alignment files will be produced
OUTDIR=BAM_strand
# MRG: This is the location where temporary
# alignment files will be created for merging.
MRG="tmp/merge"
# REFGENOME: This is the location of a .genome file
# used for bedtools.
REFGENOME=reference/CAC.genome
# CORES: This is the number of processor cores to
# be used for parallelization.
CORES=6
# EXPRDIR: This is where the bed and bedgraph coverage
# files will be produced
EXPRDIR=Expression
# CIRC: This is the location where all of the circos
# plots and configuration files will be produced.
CIRC=/home/mrals/ETP/circos
# BASE: This is the home directory of the script.
base=/home/mrals/ETP
export base
# PICARD: This is the location of picard jar files.
PICARD=/usr/local/picard-tools-1.67
# ANNOTATION: This is the location of a bed format
# CDS annotation. May be deprecated.
ANNOTATION=reference/CAC.bed
FILES=`/usr/bin/ls $INDIR/*.3.bam`

#################    P A R T   1   ###################
for file in $FILES
do
    echo $file
    f=${file##*/};f=${f%*.*.*}
    #rm $MRG/*
    ############################### P A I R E D ########################
    ########      P L U S
    # Discordant alignment type 1: dis1
    # this first type of alignment refers to a discordant type alignment
    # to the positive strand, where the first read aligns to the + strand
    # and the second read does not align to either strand (discordant)

    #samtools view -hb -f 73 $file | samtools view -hb -F 20 - > $MRG/pair.plus.for.dis1
    #samtools view -hb -f 133 $file | samtools view -hb -F 40 - > $MRG/pair.plus.rev.dis1
    # Discordant alignment type 2: dis 2
    # The second type of discordant alignment is where the first read does not
    # align to either strand and the second read aligns to the minus strand.
    #samtools view -hb -f 101 $file | samtools view -hb -F 8 - > $MRG/pair.plus.for.dis2
    #samtools view -hb -f 153 $file | samtools view -hb -F 36 - > $MRG/pair.plus.rev.dis2
    # Concordant alignment
    # This type of alignment is where the first read aligns to the forward strand
    # and the second read aligns to ther second strand, indicating that the sequenced
    # fragment originates from the forward (+) strand.
    #samtools view -hb -f 97 $file | samtools view -hb -F 28 - > $MRG/pair.plus.for.con
    #samtools view -hb -f 145 $file | samtools view -hb -F 44 - > $MRG/pair.plus.rev.con
    
    ########      M I N U S
    # Discordant alignment type 1: dis1
    # this first type of alignment refers to a discordant type alignment
    # to the negative strand, where the first read aligns to the - strand
    # and the second read does not align to either strand (discordant)

    #samtools view -hb -f 89 $file | samtools view -hb -F 36 - > $MRG/pair.minus.for.dis1
    #samtools view -hb -f 165 $file | samtools view -hb -F 24 - > $MRG/pair.minus.rev.dis1
    # Discordant alignment type 2: dis 2
    # The second type of discordant alignment is where the first read does not
    # align to either strand and the second read aligns to the + strand.
    #samtools view -hb -f 69 $file | samtools view -hb -F 40 - > $MRG/pair.minus.for.dis2
    #samtools view -hb -f 137 $file | samtools view -hb -F 52 - > $MRG/pair.minus.rev.dis2
    # Concordant alignment
    # This type of alignment is where the first read aligns to the - strand
    # and the second read aligns to the + strand, indicating that the sequenced
    # fragment originates from the reverse (-) strand.
    #samtools view -hb -f 81 $file | samtools view -hb -F 44 - > $MRG/pair.minus.for.con
    #samtools view -hb -f 161 $file | samtools view -hb -F 28 - > $MRG/pair.minus.rev.con

    ############# U N P A I R E D #########
    #samtools view -hb -F 21 $file > $MRG/un.plus
    #samtools view -hb -F 5 $file | samtools view -hb -f 16 - > $MRG/un.minus
    # Merge - sort
    plus=`/usr/bin/ls $MRG/*.plus*`
    minus=`/usr/bin/ls $MRG/*.minus*`
    #samtools merge -f tmp/plus.bam $plus 
    #samtools sort tmp/plus.bam $OUTDIR/$f.plus
    #samtools merge -f tmp/minus.bam $minus
    #samtools sort tmp/minus.bam $OUTDIR/$f.minus
    #rm tmp/plus.bam tmp/minus.bam
    #samtools index $OUTDIR/$f.plus.bam $OUTDIR/$f.plus.bam.bai TMP_DIR=tmp
    #samtools index $OUTDIR/$f.minus.bam $OUTDIR/$f.minus.bam.bai TMP_DIR=tmp

    #      P R O C E S S   B A M   F O R  C O V E R A G E


    #          B E D G R A P H

    #java -jar $PICARD/MarkDuplicates.jar INPUT=$OUTDIR/$f.plus.bam OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | samtools view -b -h -q 20 - | bedtools genomecov -bg -ibam /dev/stdin -g $REFGENOME > $EXPRDIR/$f+.bedgraph
    #java -jar $PICARD/MarkDuplicates.jar INPUT=$OUTDIR/$f.minus.bam OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | samtools view -b -h -q 20 - | bedtools genomecov -bg -ibam /dev/stdin -g $REFGENOME > $EXPRDIR/$f-.bedgraph

    #          C O V E R A G E

    #java -jar $PICARD/MarkDuplicates.jar INPUT=$OUTDIR/$f.plus.bam OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | samtools view -b -h -q 20 - | bedtools genomecov -d -ibam /dev/stdin -g $REFGENOME > $EXPRDIR/$f+.cov
    #java -jar $PICARD/MarkDuplicates.jar INPUT=$OUTDIR/$f.minus.bam OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | samtools view -b -h -q 20 - | bedtools genomecov -d -ibam /dev/stdin -g $REFGENOME > $EXPRDIR/$f-.cov

    #          G E N E   C O V E R A G E


done

#################    P A R T   2   ###################
declare -a array=()
let x=0
for file in $FILES
do
    ff=${file##*/};f=${ff%*.*.*}
    # This produces a bedpe file that omits unpaired reads and duplicates.
    #java -jar $PICARD/MarkDuplicates.jar INPUT=$file OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | samtools view -h -q 20 - | java -jar $PICARD/SortSam.jar INPUT=/dev/stdin OUTPUT=/dev/stdout SORT_ORDER=queryname | bamToBed -bedpe -mate1 > $EXPRDIR/$f.bedpe
    # This option includes duplicates
    #samtools sort -on $file - | bamToBed -bedpe -mate1 > $EXPRDIR/$f.bedpe
    # This produces a bed file that omits paired reads and duplicates.
    #java -jar $PICARD/MarkDuplicates.jar INPUT=$file OUTPUT=/dev/stdout METRICSX_FILE=/dev/null REMOVE_DUPLICATES=true | samtools view -h -q 20 -F 1 - | java -jar $PICARD/SortSam.jar INPUT=/dev/stdin OUTPUT=/dev/stdout SORT_ORDER=queryname | bamToBed > $EXPRDIR/$f.bedunpaired
    #cat $EXPRDIR/$f.bedpe | ./bedscript.rb | sort -k 1,1 -k2,2n > $EXPRDIR/$f.bedtmp
    #cat $EXPRDIR/$f.bedtmp $EXPRDIR/$f.bedunpaired | sort -k 1,1 -k2,2n > $EXPRDIR/$f.bed
    #rm $EXPRDIR/$f.bedtmp $EXPRDIR/$f.bedunpaired
    #bedtools genomecov -bg -strand + -i $EXPRDIR/$f.bed -g $REFGENOME > $EXPRDIR/$f+.bedgraph
    #bedtools genomecov -bg -strand - -i $EXPRDIR/$f.bed -g $REFGENOME > $EXPRDIR/$f-.bedgraph
    
    sed -i "s/gi|15893298|ref|NC_003030.1|/C/g" $EXPRDIR/$f+.bedgraph
    sed -i 's/gi|15004705|ref|NC_001988.2|/P/g' $EXPRDIR/$f+.bedgraph
    sed -i "s/gi|15893298|ref|NC_003030.1|/C/g" $EXPRDIR/$f-.bedgraph
    sed -i 's/gi|15004705|ref|NC_001988.2|/P/g' $EXPRDIR/$f-.bedgraph
    cp $EXPRDIR/$f+.bedgraph $EXPRDIR/$f+log.bedgraph
    cp $EXPRDIR/$f-.bedgraph $EXPRDIR/$f-log.bedgraph
    ruby -n -i -e 'temp=$_.chomp.split; puts (temp[0..2].join("\t")+"\t#{Math.log10(temp[3].to_i)}")' $EXPRDIR/$f+log.bedgraph
    ruby -n -i -e 'temp=$_.chomp.split; puts (temp[0..2].join("\t")+"\t#{Math.log10(temp[3].to_i)}")' $EXPRDIR/$f-log.bedgraph
    OLD1="file=plot1.old"
    OLD2="file=plot2.old"
    TMP="file="
    mkdir $CIRC/$f
    cp $CIRC/*.conf $CIRC/$f/
    TEMP=$TMP$base/$EXPRDIR/$f+log.bedgraph
    sed -i "s|$OLD1|$TEMP|" $CIRC/$f/circos.conf
    TEMP=$TMP$base/$EXPRDIR/$f-log.bedgraph
    sed -i "s|$OLD2|$TEMP|" $CIRC/$f/circos.conf
    array[$x]=$CIRC/$f/
    let x=$x+1
    #coverageBed -s -hist -a $EXPRDIR/$f.bed -b $ANNOTATION > $EXPRDIR/$f.counts_histogram
done
echo 'this is the array:'

parallel -j$CORES 'cd {}; circos; cd $base' ::: ${array[@]}


#samtools mpileup -BQ0 run.sorted.bam | perl -pe '($c, $start, undef, $depth) = split;if ($c ne $lastC || $start != $lastStart+1) {print "fixedStep chrom=$c start=$start step=1 span=1\n";}$_ = $depth."\n";($lastC, $lastStart) = ($c, $start);' | gzip -c > run.wig.gz





## EOF-------------------------------------------
