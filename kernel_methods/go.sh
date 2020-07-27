#!/bin/bash -ex

python3 main.py REDDIT-BINARY walk
python3 main.py REDDIT-BINARY shortest
python3 main.py REDDIT-MULTI-5K wl
python3 main.py REDDIT-MULTI-5K shortest
python3 main.py REDDIT-MULTI-5K walk
python3 main.py IMDB-BINARY wl
python3 main.py IMDB-MULTI wl
python3 main.py IMDB-BINARY shortest
python3 main.py IMDB-MULTI shortest
python3 main.py IMDB-BINARY walk
python3 main.py IMDB-MULTI walk
python3 main.py REDDIT-MULTI-5K shortest


