ROOTDIR= ../../
include $(ROOTDIR)/Makefile.base

CFLAGS= -std=c++11 -fPIC
INCLUDE= -I$(INCLUDEDIR) -I$(UTILITYDIR) 

LIBS= -lstdc++ -lm -lpthread
OBJFILES= $(UTILITYDIR)/*.o $(FILETYPESDIR)/*.o $(MATRIXDIR)/*.o 
LMACROS=


ifneq (,$(findstring withMatio,$(MACROS)))
	INCLUDE+= -I$(MATIO)/include
	LIBS+= -L$(MATIO)/lib -lmatio
	LMACROS+=-D_withMatio
endif


filetypes.o: MatrixIF.hpp *.hpp *.cpp  
	$(CC) -c filetypes.cpp $(CFLAGS) $(INCLUDE) $(LMACROS) 


objects: filetypes.o 

all: objects

clean: 
	@rm -f *.o 

anew: clean all 




# ---- DEPRECATED ----------------------------------------------------------------------------------------


MatrixIF_ASCII.o: MatrixIF.hpp MatrixIF_ASCII.hpp MatrixIF_ASCII.cpp 
	$(CC) $(LMACROS) -c MatrixIF_ASCII.cpp $(CFLAGS) $(INCLUDE)

MatrixOF_ASCII.o: MatrixOF.hpp MatrixOF_ASCII.hpp MatrixOF_ASCII.cpp 
	$(CC) $(LMACROS) -c MatrixOF_ASCII.cpp $(CFLAGS) $(INCLUDE)


MatrixIF_Matlab.o: MatrixIF.hpp MatrixIF_Matlab.hpp MatrixIF_Matlab.cpp 
	$(CC) $(LMACROS) -c MatrixIF_Matlab.cpp $(CFLAGS) $(INCLUDE)

MatrixOF_Matlab.o: MatrixOF.hpp MatrixOF_Matlab.hpp MatrixOF_Matlab.cpp 
	$(CC) $(LMACROS) -c MatrixOF_Matlab.cpp $(CFLAGS) $(INCLUDE)


MatrixIF_Boeing.o: MatrixIF.hpp MatrixIF_Boeing.hpp MatrixIF_Boeing.cpp 
	$(CC) $(LMACROS) -c MatrixIF_Boeing.cpp $(CFLAGS) $(INCLUDE)

MatrixOF_Boeing.o: MatrixOF.hpp MatrixOF_Boeing.hpp MatrixOF_Boeing.cpp 
	$(CC) $(LMACROS) -c MatrixOF_Boeing.cpp $(CFLAGS) $(INCLUDE)

#ifneq (,$(findstring withMatio,$(MACROS)))
#	$(MAKE) MatrixIF_Matlab.o 
#	$(MAKE) MatrixOF_Matlab.o 
#endif 


#objects: MatrixIF_ASCII.o MatrixOF_ASCII.o MatrixIF_Boeing.o MatrixOF_Boeing.o 

#unity: filetypes.o 



