ROOTDIR=..
include $(ROOTDIR)/Makefile.base

CFLAGS+= -fPIC#-std=c++11
INCLUDE= -I$(INCLUDEDIR) -I$(UTILITYDIR) -I$(FILETYPESDIR) #-I$(MATRIXDIR) 
OBJFILES= $(UTILITYDIR)/*.o $(FILETYPESDIR)/*.o ./*.o  
LIBS= -lstdc++ -lm 
LMACROS=

ifneq (,$(findstring withEigen,$(MACROS)))
	INCLUDE+= -I$(EIGENDIR)
	LMACROS+=-D_withEigen
else ifneq (,$(findstring withLapack,$(MACROS)))
	INCLUDE+= -I$(LAPACKDIR)/LAPACKE/include/ -I$(CBLASDIR)/include
	OBJFILES+= $(LAPACKDIR)/liblapacke.a $(LAPACKDIR)/liblapack.a $(BLASDIR)/blas_LINUX.a $(CBLASDIR)/lib/cblas_LINUX.a 
	LIBS+= -L$(FORTRANDIR)/lib/ -lgfortran
	LMACROS+=-D_withLapack
endif


matrices.o: *.hpp *.cpp 
	$(CC) -c matrices.cpp $(CFLAGS) $(INCLUDE) $(LMACROS)


objects: matrices.o 

tests: 

all: objects 

clean: 
	@rm -f *.o 

anew: clean all 
