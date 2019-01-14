#!/bin/bash -ex

DS=$1
# run preprocess
#python preprocess.py gitgraph-proc
#python preprocess.py IMDB-BINARY
#python preprocess.py IMDB-MULTI
#python preprocess.py COLLAB
#python preprocess.py DD
#python preprocess.py REDDIT-BINARY
#python preprocess.py REDDIT-MULTI-5K

~/ENV2/bin/python main.py -c ../data/$DS -l ../data/$DS.Labels -d 1024 --wlk_h 3 -e 1000 -lr 0.01
~/ENV2/bin/python main.py -c ../data/$DS -l ../data/$DS.Labels -d 512 --wlk_h 3 -e 1000 -lr 0.01
