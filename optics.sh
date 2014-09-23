#!/bin/bash
#PBS -N optics
#PBS -r n
#PBS -V
#PBS -l nodes=1:ppn=3
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
CORES=3

# DATA, MINPTS, REACH (minimum reachability distance), MAXIMA (maxima ratio),
# and DISTS are options for the optics algorithm, described @ http://github.com/MatthewRalston/OPTICS-Automatic-Clustering
declare -a DATA=('raw' 'scaled' 'correlation')
FILES=(`echo ${DATA[@]} | ruby -e 'l=gets.chomp; puts(l.split.map{|x|x+"/"+x+".csv"}.join(" "))'`)
declare -a MINPTS=({2..9})
declare -a REACH=(0.0{0..9}{1..9})
declare -a MAXIMA=(0.{5..9}{0..9})
declare -a DISTS=('braycurtis' 'canberra' 'chebyshev' 'cityblock' 'correlation' 'cosine' 'euclidean' 'mahalanobis' 'minkowski' 'seuclidean' 'sqeuclidean');




parallel -j$CORES  --tmpdir /home/mrals/Final/clustering/output --files 'optics' ::: ${DATA[@]} ::: ${MINPTS[@]} ::: ${REACH[@]} ::: ${MAXIMA[@]} ::: ${DISTS[@]}

