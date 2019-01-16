ROOTDIR= ../
include $(ROOTDIR)/Makefile.base

CFLAGS= -std=c++11 -I $(INCLUDEDIR) -I $(UTILITYDIR) -I $(FILETYPESDIR) -I $(MATRIXDIR) -fPIC
LIBS= -lstdc++ -lm -lpthread

ALLOBJECTS= $(UTILITYDIR)/*.o $(MATRIXDIR)/*.o  *.o    #$(BLOCKEDDIR)/*.o 

FLGinstance.o: FLGinstance.hpp FLGinstance.cpp
	$(CC) -c FLGinstance.cpp $(CFLAGS)

FLGkernel.o: FLGkernel.hpp FLGkernel.cpp FLGinstance.hpp
	$(CC) -c FLGkernel.cpp $(CFLAGS)

MLGgraph.o: MLGgraph.hpp MLGgraph.cpp FLGinstance.hpp FLGkernel.hpp Linearizer.hpp
	$(CC) -c MLGgraph.cpp $(CFLAGS)

MLGdataset.o: MLGdataset.hpp MLGdataset.cpp FLGinstance.hpp FLGkernel.hpp Linearizer.hpp MLGgraph.hpp
	$(CC) -c MLGdataset.cpp $(CFLAGS)

runMLG: runMLG.cpp params.hpp MLGgraph.o MLGdataset.o FLGinstance.o
	$(CC) -o runMLG runMLG.cpp $(CFLAGS) $(ALLOBJECTS) $(LIBS)

all: FLGinstance.o FLGkernel.o MLGgraph.o MLGdataset.o runMLG

clean:
	@rm -f runMLG *.o

anew: clean all
