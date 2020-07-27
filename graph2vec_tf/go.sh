#!/bin/bash -ex

# run preprocess
#python preprocess.py IMDB-BINARY
#python preprocess.py IMDB-MULTI
#python preprocess.py COLLAB
#python preprocess.py DD
#python preprocess.py REDDIT-BINARY
#python preprocess.py REDDIT-MULTI-5K

for i in 1 2 3 4 5 
do
  for DS in 'MUTAG' 'PTC_MR' 'PROTEINS_full' 'IMDB-BINARY' 'IMDB-MULTI' 'REDDIT-BINARY' 'REDDIT-MULTI-5K'
  do
  python3 preprocess.py $DS
  main.py -c ../data/$DS -l ../data/$DS.Labels -d 512 --wlk_h 3 -e 1000 -lr 0.001
  main.py -c ../data/$DS -l ../data/$DS.Labels -d 512 --wlk_h 3 -e 1000 -lr 0.01
  main.py -c ../data/$DS -l ../data/$DS.Labels -d 512 --wlk_h 3 -e 1000 -lr 0.1
  main.py -c ../data/$DS -l ../data/$DS.Labels -d 512 --wlk_h 3 -e 1000 -lr 0.5
  done
done
