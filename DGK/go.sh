#!/bin/bash -ex
DS=gitgraph-proc
source activate my-rdkit-env2
#python deep_kernel.py 512 3 3 $DS 5 1 1 7 100
#python deep_kernel.py 512 3 3 $DS 5 1 1 7 100
#python deep_kernel.py 512 3 3 $DS 5 1 1 7 100
#python deep_kernel.py 512 3 3 $DS 5 1 1 7 100
#python deep_kernel.py 512 3 3 $DS 5 1 1 7 100

#python2 deep_kernel.py 512 3 1 $DS 5 1 1 7 100
#python2 deep_kernel.py 512 3 1 $DS 5 1 1 7 100
#python2 deep_kernel.py 512 3 1 $DS 5 1 1 7 100
#python2 deep_kernel.py 512 3 1 $DS 5 1 1 7 100
#python2 deep_kernel.py 512 3 1 $DS 5 1 1 7 100

python2 deep_kernel.py 512 3 2 $DS 5 1 1 7 100
python2 deep_kernel.py 512 3 2 $DS 5 1 1 7 100
python2 deep_kernel.py 512 3 2 $DS 5 1 1 7 100
python2 deep_kernel.py 512 3 2 $DS 5 1 1 7 100
python2 deep_kernel.py 512 3 2 $DS 5 1 1 7 100
