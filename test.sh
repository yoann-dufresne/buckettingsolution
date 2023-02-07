#!/bin/bash

DATA=data/sarscov2.fa

K=31
M=11

set -e
make

OUT="${DATA%.*}_k${K}_m$M"
rm -rf $OUT
mkdir $OUT

echo "./bucketting $DATA $OUT $K $M"
/usr/bin/time -v ./bucketting $DATA $OUT $K $M

DATAN="$DATA.noN.fa"
cp $DATA $DATAN
sed -i 's/N//g' $DATAN
python3 verif.py $DATAN $OUT $K

rm $DATAN
