#!/bin/bash
#PBS -N visualization
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -d /home/mrals/ETP

#set -e

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
PICARD=/usr/local/picard-tools-1.67

#########################
# S C R I P T
#########################

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
    echo 'mergensort'
    plus=`/usr/bin/ls $MRG/*.plus*`
    minus=`/usr/bin/ls $MRG/*.minus*`
    #samtools merge -f tmp/plus.bam $plus 
    #samtools sort tmp/plus.bam $OUTDIR/$f.plus
    #samtools merge -f tmp/minus.bam $minus
    #samtools sort tmp/minus.bam $OUTDIR/$f.minus
    #rm tmp/plus.bam tmp/minus.bam
    #samtools index $OUTDIR/$f.plus.bam $OUTDIR/$f.plus.bam.bai TMP_DIR=tmp
    #samtools index $OUTDIR/$f.minus.bam $OUTDIR/$f.minus.bam.bai TMP_DIR=tmp
done

#################    P A R T   2   ###################
declare -a array=()
let x=0
for file in $FILES
do
    f=${file##*/};f=${f%*.*.*}
    # This produces a bedpe file that omits unpaired reads and duplicates.
    java -jar $PICARD/MarkDuplicates.jar INPUT=$file OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | samtools view -h -q 20 - | java -jar $PICARD/SortSam.jar INPUT=/dev/stdin OUTPUT=/dev/stdout SORT_ORDER=queryname | bamToBed -bedpe -mate1 > $EXPRDIR/$f.bedpe
    # This option includes duplicates
    #samtools sort -on $file - | bamToBed -bedpe -mate1 > $EXPRDIR/$f.bedpe
    # This produces a bed file that omits paired reads and duplicates.
    java -jar $PICARD/MarkDuplicates.jar INPUT=$file OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | samtools view -h -q 20 -F 1 - | java -jar $PICARD/SortSam.jar INPUT=/dev/stdin OUTPUT=/dev/stdout SORT_ORDER=queryname | bamToBed > $EXPRDIR/$f.bedunpaired
    cat $EXPRDIR/$f.bedpe | ./bedscript.rb | sort -k 1,1 -k2,2n > $EXPRDIR/$f.bedtmp
    cat $EXPRDIR/$f.bedtmp $EXPRDIR/$f.bedunpaired | sort -k 1,1 -k2,2n > $EXPRDIR/$f.bed
    rm $EXPRDIR/$f.bedtmp $EXPRDIR/$f.bedunpaired
    bedtools genomecov -bg -strand + -i $EXPRDIR/$f.bed -g $REFGENOME > $EXPRDIR/$f+.cov
    bedtools genomecov -bg -strand - -i $EXPRDIR/$f.bed -g $REFGENOME > $EXPRDIR/$f-.cov
    sed -i "s/gi|15893298|ref|NC_003030.1|/C/g" $EXPRDIR/$f+.cov
    sed -i 's/gi|15004705|ref|NC_001988.2|/P/g' $EXPRDIR/$f+.cov
    sed -i "s/gi|15893298|ref|NC_003030.1|/C/g" $EXPRDIR/$f-.cov
    sed -i 's/gi|15004705|ref|NC_001988.2|/P/g' $EXPRDIR/$f-.cov
    #OLD1="file=plot1.old"
    #OLD2="file=plot2.old"
    #TMP="file="
    #mkdir $CIRC/$f
    #cp $CIRC/*.conf $CIRC/$f/
    #TEMP=$TMP$base/$EXPRDIR/$f+.cov
    #sed -i "s|$OLD1|$TEMP|" $CIRC/$f/circos.conf
    #TEMP=$TMP$base/$EXPRDIR/$f-.cov
    #sed -i "s|$OLD2|$TEMP|" $CIRC/$f/circos.conf
    #array[$x]=$base/$CIRC/$f/
    #let x=$x+1

done

#parallel -j$CORES 'cd {}; circos; cd $base' ::: $array







## EOF-------------------------------------------
