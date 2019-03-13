#!/bin/bash -ex
for i in 1 2 3 4 5
do
  for DS in 'IMDB-BINARY' 'IMDB-MULTI'
  do
   python3 src/main.py --input ../data/$DS --preprocessed-input preprocessed_dataset/$DS --d 512 --property n
   python3 src/main.py --input ../data/$DS --preprocessed-input preprocessed_dataset/$DS --d 512 --property s 
  done
done
