from MLGK import *
import numpy as np
import pdb
import os.path as path
import sys

fname =  path.abspath(path.join(__file__ ,"../../../data/MUTAG.txt"))
save_name = "m.txt"
num_graphs = 188
eta = 0.01
gamma= 0.001
grow = True
radius = 2
levels = 2

m = MLGdataset(fname, eta, gamma, grow)
m.computeGram(levels, radius)
m.saveGram(save_name)

# fill a numpy matrix. Note: need to call compute gram
# before calling fillGram
gram = np.zeros((num_graphs, num_graphs))
m.fillGram(gram)
pdb.set_trace()
