#!/bin/bash

INPUT=$1
CONSENSUS=$2
OUTPUT=$3

ID=${INPUT%den_input}

DIR=/projects/ps-gleesonlab8/User/hiyoothere/NTD/10.Poisson
 
module load bedtools
python $DIR/locus_to_bed.py $INPUT  $ID'.bed'
sort -k1,1V -k2,2g $ID'.bed' > $ID'sorted.bed'
bedtools intersect -a  $ID'sorted.bed' -b $CONSENSUS >  $ID'sorted.consensus.bed'
python $DIR/bed_to_locus.py $ID'sorted.consensus.bed'  $INPUT  $OUTPUT

rm $ID'.bed'
rm $ID'sorted.bed'
rm $ID'sorted.consensus.bed'
