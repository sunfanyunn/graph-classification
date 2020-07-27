#!/bin/bash -ex

# Fill in the name of the dataset
DS=

# Run multiple trials 
python2 deep_kernel.py 512 3 2 $DS 5 1 1 7 100
python2 deep_kernel.py 512 3 2 $DS 5 1 1 7 100
python2 deep_kernel.py 512 3 2 $DS 5 1 1 7 100
python2 deep_kernel.py 512 3 2 $DS 5 1 1 7 100
python2 deep_kernel.py 512 3 2 $DS 5 1 1 7 100
