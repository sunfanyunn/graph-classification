#!/bin/bash -ex
# Run the MLG kernel on the MUTAG dataset with parameters:
# radius = 1
# levels = 2
# eta = 0.1
# gamma = 0.01
# num threads = 32
# grow = 1 # if you want the subgraphs to double in size at each level, set this equal to 0

# Replace MUTAG with the dataset name of your choice(PTC/PROTEINS/NCI1/NCI109).
BASE=`pwd`
dset=$1
data=$BASE/../data/$dset.txt
feats=$BASE/..//data/$dset\_nodelabels.txt
save=$BASE//data/results/output.txt
mkdir -p $BASE/data/results/

~/ENV/bin/python3 preprocess.py $dset

for r in 1 2 3 4
do 
  for l in 1 2 3 4
  do
    for g in 0.01 0.1 1
    do
      for e in 0.01 0.1 1
      do

  cd MLGkernel
  ./runMLG -d $data -f $feats -s $save -r $r -l $l -e $e -g $g -t 32 -m 1
  cd ../
  ~/ENV/bin/python3 evaluate_embedding.py $dset >> $dset.log
      done
    done
  done
done
