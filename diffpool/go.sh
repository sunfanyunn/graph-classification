#!/bin/bash -ex
# grid search for gitgraph-proc
for method in 'base-set2set' 'base'
  do
  for gc in 4 8 16 2
  do
    python -m train --datadir=../data --bmname=gitgraph-proc --cuda=0 --max-nodes 1000 --epochs=100 --num-classes=3 --output-dim 512 --lr 0.001 --num-gc-layers $gc --method $method 
  done
done

# DD
#python -m train --datadir=data --bmname=DD --cuda=0 --max-nodes=500 --epochs=1000

# ENZYMES
#python -m train --datadir=data --bmname=ENZYMES --cuda=3 --max-nodes=100 --num-classes=6

