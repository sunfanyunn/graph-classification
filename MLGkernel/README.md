# The Multiscale Laplacian Graph Kernel

This is a C++ implementation of the Multiscale Graph Laplacian kernel as described in:
 
R. Kondor, H. Pan, [The Multiscale Graph Laplacian](https://arxiv.org/abs/1603.06186) (2016)

## Requirements
* C++11
* [Eigen](http://eigen.tuxfamily.org/index.php)

## Installation/Setup
Change the EIGENDIR variable Makefile.options to the path to your installation of
the Eigen library. Run the following command to create the runMLG executable in the MLGkernel directory.
```bash
make all
```

## Run the demo
After compiling, you can run the sample.sh script, which will run the MLGkernel
on the sample datasets in the data directory and save the resulting kernel matrices in data/results.
```bash
sh sample.sh
```

## Content
To understand the code, the main files worth skimming through are in the MLGkernel directory:
MLGdataset, Linearizer, FLGinstance, and FLGkernel.
See MLGkernel/runMLG.cpp for an example of how to load data/feature files and
how to compute/save the MLG kernel matrix.

## Data Format
In order to use your own data, you have to provide the data via
* a text file in the following format: the first line is the number of graphs in the dataset
The next line is the number of vertices in first graph, the next {size of first graph} lines
is the space delimited adjacency matrix of the first graph. And so on for the rest of the graphs.
See data/MUTAG.txt for an example of this format.
* a text file containing a discrete node label of each graph for each graph in the dataset.
The first line is the number of graphs in the dataset. The second line is the size of the first
graph in the dataset. The next {size of first graph} lines denote the node labels of the vertices
of the first graph. And so on for the rest of the graphs.
See data/MUTAG_nodelabels.txt for an example of this data format.

Note that if no feature file is given for the runMLG
executable, it will default to using vertex degrees as node labels.

## Contact
If you have any questions/comments/concerns, please contact hopan[at]uchicago.edu. In particular,
it would be fantastic if you could report any bugs!
