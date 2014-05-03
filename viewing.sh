#!/bin/bash
#PBS -N visualization
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -d /home/mrals/ETP

set -e

. ~/.bash_profile
#General
PATH=$PATH:/home/mrals/pckges/Vienna/bin:/home/mrals/home/bin/:/home/mrals/bin/
#NGS
PATH=$PATH:/home/mrals/pckges/sratoolkit.2.3.4-2-centos_linux64/bin/:/usr/local/bowtie-0.12.7/:/usr/local/bowtie2-2.1.0/:/usr/local/bwa-0.7.4/:/usr/local/cufflinks-2.0.2/:/usr/local/FastQC/:/usr/local/samtools-0.1.18/:/usr/local/tophat-2.0.4/:/home/mrals/pckges/seqtk-master/

export PATH
[[ -s "$HOME/.rvm/scripts/rvm" ]] && source "$HOME/.rvm/scripts/rvm"
rvm use 2.0.0

#########################
#  D E S C R I P T I O N
#########################
#  AUTHOR: MATT RALSTON
#  DATE: 1/22/14
#
# This script contains code that is divided into two parts.
#
#################    P A R T   1   ###################
# The first part is designed to separate processed BAM files by strand.
# Condordant, discordant, and unpaired strand-speicifc alignments will be
# separated into distinct files, for visualization in a genomic viewer.

################     P A R T   2   ###################
# The second part is designed to visualize the coverage of processed BAM files
# by strand in a circos plot.


#########################
# ARGUMENTS
#########################
# Optional for RNAseq: RefFlat file
# REFFLAT=reference/H.pylori.refFlat
INDIR=SAM_processed
OUTDIR=BAM_strand
FILES=`/usr/bin/ls $INDIR/*.3.bam`
MRG="tmp/merge"
REFGENOME=reference/CAC.genome
CORES=6
EXPRDIR=Expression
CIRC=/home/mrals/ETP/circos
base=/home/mrals/ETP

#########################
# S C R I P T
#########################

#################    P A R T   1   ###################
#for file in $FILES
#do
#    rm $MRG/*
    ############################### P A I R E D ########################
    ########      P L U S
    # Discordant alignment type 1: dis1
    # this first type of alignment refers to a discordant type alignment
    # to the positive strand, where the first read aligns to the + strand
    # and the second read does not align to either strand (discordant)
#    samtools view -f 73 $file - | samtools view -F 20 - > $MRG/pair.plus.for.dis1
#    samtools view -f 133 $file - | samtools view -F 40 - > $MRG/pair.plus.rev.dis1
    # Discordant alignment type 2: dis 2
    # The second type of discordant alignment is where the first read does not
    # align to either strand and the second read aligns to the minus strand.
#    samtools view -f 101 $file - | samtools view -F 8 - > $MRG/pair.plus.for.dis2
#    samtools view -f 153 $file - | samtools view -F 36 - > $MRG/pair.plus.rev.dis2
    # Concordant alignment
    # This type of alignment is where the first read aligns to the forward strand
    # and the second read aligns to ther second strand, indicating that the sequenced
    # fragment originates from the forward (+) strand.
#    samtools view -f 97 $file - | samtools view -F 28 - > $MRG/pair.plus.for.con
#    samtools view -f 145 $file - | samtools view -F 44 - > $MRG/pair.plus.rev.con

    ########      M I N U S
    # Discordant alignment type 1: dis1
    # this first type of alignment refers to a discordant type alignment
    # to the negative strand, where the first read aligns to the - strand
    # and the second read does not align to either strand (discordant)
#    samtools view -f 89 $file - | samtools view -F 36 - > $MRG/pair.minus.for.dis1
#    samtools view -f 165 $file - | samtools view -F 24 - > $MRG/pair.minus.rev.dis1
    # Discordant alignment type 2: dis 2
    # The second type of discordant alignment is where the first read does not
    # align to either strand and the second read aligns to the + strand.
#    samtools view -f 69 $file - | samtools view -F 40 - > $MRG/pair.minus.for.dis2
#    samtools view -f 137 $file - | samtools view -F 52 - > $MRG/pair.minus.rev.dis2
    # Concordant alignment
    # This type of alignment is where the first read aligns to the - strand
    # and the second read aligns to the + strand, indicating that the sequenced
    # fragment originates from the reverse (-) strand.
#    samtools view -f 81 $file - | samtools view -F 44 - > $MRG/pair.minus.for.con
#    samtools view -f 161 $file - | samtools view -F 28 - > $MRG/pair.minus.rev.con

    ############# U N P A I R E D #########
#    samtools view -F 21 $file $MRG/un.plus
#    samtools view -F 5 $file - | samtools view -f 16 - > $MRG/un.minus
    # Merge - sort
#    plus=`/usr/bin/ls $MRG/*.plus*`
#    minus=`/usr/bin/ls $MRG/*.minus*`
#    samtools merge - $plus | samtools sort - $OUTDIR/${file%*.*.*}.plus.bam
#    samtools merge - $minus | samtools sort - $OUTDIR/${file%*.*.*}.minus.bam
#done

#################    P A R T   2   ###################
declare -a array=()
let x=0
for file in $FILES
do
    f=${file##*/};f=${f%*.*.*}
    #samtools sort -on $file - | bamToBed -bedpe -mate1 > $EXPRDIR/$f.bedpe 
    #cat $EXPRDIR/$f.bedpe | ./bedscript.rb | sort -k 1,1 > $EXPRDIR/$f.bed
    #bedtools genomecov -bga -strand + -i $EXPRDIR/$f.bed -g $REFGENOME > $EXPRDIR/$f+.cov
    #bedtools genomecov -bga -strand - -i $EXPRDIR/$f.bed -g $REFGENOME > $EXPRDIR/$f-.cov
    grep -rl "gi|15893298|ref|NC_003030.1|" $base/$EXPRDIR/$f+.cov | xargs sed -i 's/gi|15893298|ref|NC_003030.1|/C/g'
    grep -rl "gi|15004705|ref|NC_001988.2|" $base/$EXPRDIR/$f+.cov | xargs sed -i 's/gi|15004705|ref|NC_001988.2|/P/g'
    grep -rl "gi|15893298|ref|NC_003030.1|" $base/$EXPRDIR/$f-.cov | xargs sed -i 's/gi|15893298|ref|NC_003030.1|/C/g'
    grep -rl "gi|15004705|ref|NC_001988.2|" $base/$EXPRDIR/$f-.cov | xargs sed -i 's/gi|15004705|ref|NC_001988.2|/P/g'
    # The file is then prepared for use in creating circos plots
    OLD1="file=plot1.old"
    OLD2="file=plot2.old"
    TMP="file="
    cp $CIRC/*.conf $CIRC/$f/
    cd $CIRC/$f
    TEMP=$TMP$EXPRDIR/$f+.cov
    grep -rl $OLD1 circos.conf | xargs sed -i 's/'$OLD1'/'$TEMP'/g'
    TEMP=$TMP$EXPRDIR/$f-.cov
    grep -rl $OLD2 circos.conf | xargs sed -i 's/'$OLD2'/'$TEMP'/g'
    array[$x]=$base/$CIRC/$f/
    let x=$x+1
    cd $base

done

#parallel -j$CORES 'cd {}; circos; cd $base' ::: $array







## EOF-------------------------------------------
