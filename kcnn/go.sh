#!/bin/bash -ex

# TEST is the virtualenv 
./TEST/bin/python3 main.py --d 512 --dataset $@
