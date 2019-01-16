Run 'make all' at the top level of the MLGkernel directory.
Then in this directory("MLGkernel/MLGkernel/swig/"), run make all.
This will create the _MLGK.so shared object file as well as
the MLGK.py file.

In python you should then be able to stuff like:
import MLGK
lvls = 2
radius = 2
m = MLGK.MLGdataset()
m.loadGraphs("sampledata.txt")
m.computeGram(lvls, radius)
m.saveGram("savefile.txt")


See test.py for further details.

NOTE: You will have to change some of the macros in the Makefile here(to change
the location of where the Python.h file is stored - it should probably be
somewhere like /usr/include/python3.x/. Also, in the top level directory inside
Makefile.options, you need to specify where the eigen library is located, as well
as the C++ compiler(IE: clang/g++/etc).
