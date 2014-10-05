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
    mrg=$MRG/$f
    #mkdir $mrg
############################### P A I R E D ########################
    ########      P L U S
    # Discordant alignment type 1: dis1
    # this first type of alignment refers to a discordant type alignment
    # to the positive strand, where the first read aligns to the + strand
    # and the second read does not align to either strand (discordant)
    #samtools view -hu -f 73 $file | samtools view -hb -F 20 -o $mrg/pair.plus.for.dis1 -
    #samtools view -hu -f 133 $file | samtools view -hb -F 40 -o $mrg/pair.plus.rev.dis1 -
    # Discordant alignment type 2: dis 2
    # The second type of discordant alignment is where the first read does not
    # align to either strand and the second read aligns to the minus strand.
    #samtools view -hu -f 101 $file | samtools view -hb -F 8 -o $mrg/pair.plus.for.dis2 -
    #samtools view -hu -f 153 $file | samtools view -hb -F 36 -o $mrg/pair.plus.rev.dis2 -
    # Concordant alignment
    # This type of alignment is where the first read aligns to the forward strand
    # and the second read aligns to ther second strand, indicating that the sequenced
    # fragment originates from the forward (+) strand.
    #samtools view -hu -f 97 $file | samtools view -hb -F 28 -o $mrg/pair.plus.for.con -
    #samtools view -hu -f 145 $file | samtools view -hb -F 44 -o $mrg/pair.plus.rev.con -
    
    ########      M I N U S
    # Discordant alignment type 1: dis1
    # this first type of alignment refers to a discordant type alignment
    # to the negative strand, where the first read aligns to the - strand
    # and the second read does not align to either strand (discordant)

    #samtools view -hu -f 89 $file | samtools view -hb -F 36 -o $mrg/pair.minus.for.dis1 -
    #samtools view -hu -f 165 $file | samtools view -hb -F 24 -o $mrg/pair.minus.rev.dis1 -
    # Discordant alignment type 2: dis 2
    # The second type of discordant alignment is where the first read does not
    # align to either strand and the second read aligns to the + strand.
    #samtools view -hu -f 69 $file | samtools view -hb -F 40 -o $mrg/pair.minus.for.dis2 -
    #samtools view -hu -f 137 $file | samtools view -hb -F 52 -o $mrg/pair.minus.rev.dis2 -
    # Concordant alignment
    # This type of alignment is where the first read aligns to the - strand
    # and the second read aligns to the + strand, indicating that the sequenced
    # fragment originates from the reverse (-) strand.
    #samtools view -hu -f 81 $file | samtools view -hb -F 44 -o $mrg/pair.minus.for.con -
    #samtools view -hu -f 161 $file | samtools view -hb -F 28 -o $mrg/pair.minus.rev.con -
    ############# U N P A I R E D #########
    #samtools view -hb -F 21 -o $mrg/un.plus $file
    #samtools view -hu -F 5 $file | samtools view -hb -f 16 -o $mrg/un.minus -
    # Merge - sort
    #plus=`/usr/bin/ls $mrg/*.plus*`
    #minus=`/usr/bin/ls $mrg/*.minus*`
    #samtools merge - $plus | samtools view -h - | ruby -ne 'line=$_.chomp; line[0..2] == "HWI" ? (line=line.split("\t"); puts(([line[0].split("_")[0]]+line[1..-1]).join("\t"));) : puts(line)' | samtools view -hbS -o $mrg/plus.bam -
    #samtools sort $mrg/plus.bam $OUTDIR/$f.plus
    #samtools merge - $minus | samtools view -h - | ruby -ne 'line=$_.chomp; line[0..2] == "HWI" ? (line=line.split("\t"); puts(([line[0].split("_")[0]]+line[1..-1]).join("\t"));) : puts(line)' | samtools view -hbS -o $mrg/minus.bam -
    #samtools sort $mrg/minus.bam $OUTDIR/$f.minus
    #rm $mrg/*
    #samtools index $OUTDIR/$f.plus.bam $OUTDIR/$f.plus.bam.bai TMP_DIR=tmp
    #samtools index $OUTDIR/$f.minus.bam $OUTDIR/$f.minus.bam.bai TMP_DIR=tmp
    #      P R O C E S S   B A M   F O R  C O V E R A G E
    #          B E D G R A P H
    #java -jar $PICARD/MarkDuplicates.jar QUIET=true INPUT=$OUTDIR/$f.plus.bam OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | 
    #samtools view -bu -q 20 $OUTDIR/$f.plus.bam | bedtools genomecov -bg -ibam /dev/stdin -g $REFGENOME > $EXPRDIR/$f+.bedgraph
    #java -jar $PICARD/MarkDuplicates.jar QUIET=true INPUT=$OUTDIR/$f.minus.bam OUTPUT=/dev/stdout METRICS_FILE=/dev/null REMOVE_DUPLICATES=true | 
    #samtools view -bu -q 20 $OUTDIR/$f.minus.bam | bedtools genomecov -bg -ibam /dev/stdin -g $REFGENOME > $EXPRDIR/$f-.bedgraph
    #          C O V E R A G E
    # option to set genomecov option to -bg (areas of continuous coverage together)
    #samtools view -b -h -q 20 $OUTDIR/$f.plus.bam | bedtools genomecov -d -ibam /dev/stdin -g $REFGENOME > $EXPRDIR/$f+.cov
    #samtools view -b -h -q 20 $OUTDIR/$f.minus.bam | bedtools genomecov -d -ibam /dev/stdin -g $REFGENOME > $EXPRDIR/$f-.cov
    #sed -i "s/NC_003030.1/C/g" $EXPRDIR/$f+.bedgraph
    #sed -i "s/NC_001988.2/P/g" $EXPRDIR/$f+.bedgraph
    #sed -i "s/NC_003030.1/C/g" $EXPRDIR/$f-.bedgraph
    #sed -i "s/NC_001988.2/P/g" $EXPRDIR/$f-.bedgraph
    cat $EXPRDIR/$f+.bedgraph | ruby -ne 'line=gets.chomp.split; puts(line[0..2].join("\t")+"\t#{Math.log10(line[3].to_i).round(4)}\tsvgcoord=#{line[1..2].join("-")},svgexp=#{line[3]}")' > $EXPRDIR/$f+.log.bedgraph
    cat $EXPRDIR/$f-.bedgraph | ruby -ne 'line=gets.chomp.split; puts(line[0..2].join("\t")+"\t#{Math.log10(line[3].to_i).round(4)}\tsvgcoord=#{line[1..2].join("-")},svgexp=#{line[3]}")' > $EXPRDIR/$f-.log.bedgraph

}

function optics {
    file="$1/$1.csv"
    minpts=$2
    reach=$3
    maxima=$4
    #clust=$5
    dist=$5
    area=$6
    n="n"
    outdir=''
    cmd="OPTICS -i $file"
    #outdir="$1/min_$2_reach_$3_maxima_$4_dist_$5"
    #cmd="OPTICS -i $file -n $minpts -r $reach -m $maxima -d $dist"
    [[ $minpts != $n ]] && cmd="$cmd -n $minpts" && outdir="${outdir}min_${minpts}"
    [[ $reach != $n ]] && cmd="$cmd -r $reach" && outdir="${outdir}reach_${reach}"
    [[ $maxima != $n ]] && cmd="$cmd -m $maxima" && outdir="${outdir}maxima_${maxima}"
    #[[ $dist != $n ]] && cmd="$cmd -c $clust" && outdir="${outdir}clust_${clust}"
    [[ $dist != $n ]] && cmd="$cmd -d $dist" && outdir="${outdir}dist_${dist}"
    [[ $area != $n ]] && cmd="$cmd -s $area" && outdir="${outdir}area_${area}"
    filename=$outdir
    outdir="${1}/${outdir}"
    mkdir $outdir
    cmd="$cmd -o $outdir"
    echo "${outdir},${1},${minpts},${reach},${maxima},${dist},${area},$(eval $cmd)"  > $outdir/summary.out
    cat $outdir/summary.out > output/$filename.out
}

function mintrin {
    # takes the condition prefix as the first argument


    prefix=$1
    # takes the fastq file input directory as the second argument
    fastq=$2
    # takes the BAM file input directory as the third argument
    bam=$3
    REFGENOME=$4
    REFPROTEOME=$5
    CDS=$6

    OUT=$(echo $prefix | ruby -e 'test=gets.chomp; puts(test.gsub(/\*/,"-"))')
    mkdir $OUT

    #zcat $fastq/${prefix}_unpaired.gz | ./../untangler.rb 1 | gzip > $OUT/u1.fq.gz
    #zcat $fastq/${prefix}_unpaired.gz | ./../untangler.rb 2 | gzip > $OUT/u2.fq.gz
    #TEX=$(ruby -e 'test=`echo $prefix`; puts("true") if test.include?("75") || test.include?("270")')

    # fastq
    #zcat $fastq/$prefix.fastq.1.gz $OUT/u1.fq.gz | ruby -ne '$_[0..3] == "@HWI" ? puts($_.chomp+"/1") : puts($_.chomp)' | gzip > $OUT/left_combined.fq.gz
    #zcat $fastq/$prefix.fastq.2.gz | ruby -ne '$_[0..3] == "@HWI" ? puts($_.chomp+"/2") : puts($_.chomp)' | gzip > $OUT/right_combined.fq.gz
    # Corrects an issue with all unpaired right reads being reverse complemented in trimfiltercheck
    #zcat $OUT/u2.fq.gz | ruby -ne '$_[0..3] == "@HWI" ? puts($_.chomp+"/2") : puts($_.chomp)' | fastx_reverse_complement -Q33 | gzip >> $OUT/right_combined.fq.gz

    #if [[ "$TEX" == true ]] ; then
#	zcat $fastq/$prefix-TEX.fastq.1.gz | ruby -ne '$_[0..3] == "@HWI" ? puts($_.chomp+"/1") : puts($_.chomp)' | gzip >> $OUT/left_combined.fq.gz
#	zcat $fastq/$prefix-TEX.fastq.2.gz | ruby -ne '$_[0..3] == "@HWI" ? puts($_.chomp+"/2") : puts($_.chomp)' | gzip >> $OUT/right_combined.fq.gz
#	zcat $fastq/$prefix-TEX_unpaired.gz | ./untangler.rb 1 | ruby -ne '$_[0..3] == "@HWI" ? puts($_.chomp+"/2") : puts($_.chomp)' | gzip >> $OUT/left_combined.fq.gz
#	zcat $fastq/$prefix-TEX_unpaired.gz | ./untangler.rb 2 | ruby -ne '$_[0..3] == "@HWI" ? puts($_.chomp+"/2") : puts($_.chomp)' | fastx_reverse_complement -Q33 | gzip >> $OUT/right_combined.fq.gz
    #fi
    #rm $OUT/u1.fq.gz $OUT/u2.fq.gz
    # Merge and sort bam files
    #samtools merge $OUT/All.merged.bam `/usr/bin/ls -d $bam/$prefix.3.bam` 
    #samtools sort -o $OUT/All.merged.bam $OUT/All.most > $OUT/All.most.bam
    # Add suffixes
    #samtools view -h $OUT/All.most.bam | ruby -ne '$_[0..2] == "HWI" ? puts($_.split[0]+"/"+$_.split("_")[1][0]+"\t"+$_.split[1..$_.split.size].join("\t")) : puts($_.chomp)' | ruby -ne '$_[0..2] == "HWI" ? (  l=$_.split("\t"); l[1].to_i%2 == 1 ? puts(l.join("\t")) : puts( ([l[0]]+[(l[1].to_i+1).to_s]+l[2..l.size]).join("\t") )  ) : puts($_.chomp)' | samtools view -hbS - > $OUT/All.bam

    # Trinity, then again if any commands fail
    #Trinity --left $OUT/left_combined.fq.gz --right $OUT/right_combined.fq.gz --genome $REFGENOME --genome_guided_max_intron 1 --jaccard_clip --genome_guided_use_bam $OUT/All.bam --JM $JM --seqType fq --output $OUT --genome_guided_CPU 4 --CPU $CORES --SS_lib_type FR >> $OUT/ref_assembly.log 2>&1
    #TEST=$(ruby -e 'test=`/usr/bin/ls $prefix/`; test.split.each {|x| puts("FAIL") if x.include?("fail") || x.include?("Fail")};')
    #[[ $TEST == "FAIL" ]] && Trinity --left $OUT/left_combined.fq.gz --right $OUT/right_combined.fq.gz --genome $REFGENOME --genome_guided_max_intron 1 --jaccard_clip --genome_guided_use_bam $OUT/All.bam --JM $JM --seqType fq --output $OUT --genome_guided_CPU 4 --CPU $CORES --SS_lib_type FR >> $OUT/ref_assembly.log 2>&1
    rm -rf $OUT/Dir_* $OUT/All* $OUT/*.fq.gz
    cd $OUT
    transrate -a Trinity-GG.fasta -r $REFPROTEOME -g $REFGENOME -x 1 -f transrate_output.csv
    cd ..
    blat $REFGENOME $OUT/Trinity-GG.fasta -maxIntron=0 $OUT/Trinity.psl
    psl2bed-best $OUT/Trinity.psl /dev/stdout | ruby -ne 'puts($_.split("\t")[0...6].join("\t"))' | bedToGenePred stdin stdout| genePredToGtf file stdin stdout | grep -v 'codon' | grep -v 'CDS' | sed 's/stdin/trinity/' | sed 's/exon\t/transcript\t/' | sed 's/|/-/' | ./assemblycurate.rb > $OUT/Trinity_filtered.gtf
    ./assemblycompare.rb $OUT/Trinity_filtered.gtf $CDS $OUT
}

export -f quality preprocess bamsplit optics mintrin

## EOF------------------------------------------
