#!/bin/bash
#PBS -N optics
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=5
#PBS -l walltime=12:00:00
#PBS -d /home/mrals/Final/
#------------------------------------------------
# Title: optics.sh
#
# Copyright 2014 Matt Ralston
#
# This script uses several different input files
# and parameter values to cluster using optics.
#------------------------------------------------


#set -e
source functions.sh
cd clustering
#------------------------------------------------
# Parameters
#------------------------------------------------
# CORES
CORES=5



# DATA, MINPTS, REACH (minimum reachability distance), MAXIMA (maxima ratio),
# and DISTS are options for the optics algorithm, described @ http://github.com/MatthewRalston/OPTICS-Automatic-Clustering
declare -a DATA=('raw' 'scaled' 'pearson' 'spearman' 'kendall')
FILES=(`echo ${DATA[@]} | ruby -e 'l=gets.chomp; puts(l.split.map{|x|x+"/"+x+".csv"}.join(" "))'`)
declare -a MINPTS=({2..8..2})
declare -a REACH=(0.0{0..8..2}{1..9..2})
declare -a MAXIMA=(0.{3..9}{0..8..2})
declare -a DISTS=('braycurtis' 'canberra' 'chebyshev' 'cityblock' 'correlation' 'hamming' 'seuclidean' 'euclidean')
declare -a AREA=(0.{6..8}{0..8..4})

#parallel -j$CORES  --tmpdir /home/mrals/Final/clustering/output --files 'optics' ::: ${DATA[@]} ::: ${MINPTS[@]} ::: ${REACH[@]} ::: ${MAXIMA[@]} ::: ${DISTS[@]}

time parallel -j$CORES  --tmpdir /home/mrals/Final/clustering/output --files 'optics' ::: ${DATA[@]} ::: ${MINPTS[@]} ::: "n" ::: ${MAXIMA[@]} ::: ${DISTS[@]} ::: ${AREA[@]}


