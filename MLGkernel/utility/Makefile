ROOTDIR= ..
include $(ROOTDIR)/Makefile.base

CFLAGS+= -I $(INCLUDEDIR) -I $(UTILITYDIR) -I $(MATRIXDIR) -fPIC # -std=c++11


Log.o: Log.hpp Log.cpp
	$(CC) -c Log.cpp $(CFLAGS)

Rstream.o: Rstream.hpp Rstream.cpp
	$(CC) -c Rstream.cpp $(CFLAGS)

Bofstream.o: Bofstream.hpp Bofstream.cpp
	$(CC) -c Bofstream.cpp $(CFLAGS)

Bifstream.o: Bifstream.hpp Bifstream.cpp
	$(CC) -c Bifstream.cpp $(CFLAGS)

Serializable.o: Serializable.hpp Serializable.cpp
	$(CC) -c Serializable.cpp $(CFLAGS)

TopkList.o: TopkList.hpp TopkList.cpp 
	$(CC) -c TopkList.cpp $(CFLAGS) 

ThreadBank.o: ThreadBank.hpp ThreadBank.cpp 
	$(CC) -c ThreadBank.cpp $(CFLAGS)

ThreadManager.o: ThreadManager.hpp ThreadManager.cpp ThreadBank.hpp 
	$(CC) -c ThreadManager.cpp $(CFLAGS)

#Graph.o: Graph.cpp Graph.hpp
#	$(CC) -c Graph.cpp $(CFLAGS) $(INCLUDE)

filetypes.o: filetypes/*.hpp filetypes/*.cpp  
	$(CC) -c filetypes/filetypes.cpp $(CFLAGS) $(INCLUDE) -I filetypes/ $(MACROS) 


objects: Log.o Rstream.o Bofstream.o Bifstream.o Serializable.o TopkList.o \
ThreadBank.o ThreadManager.o filetypes.o #Graph.o

tests:

all: objects tests 

clean: 
	@rm -f *.o 

anew: clean all 
