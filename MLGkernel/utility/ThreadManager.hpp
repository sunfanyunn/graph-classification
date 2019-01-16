/*
 -----------------------------------------------------------------------------
 
 MLGkernel is an open source implementation of the Multiscale Laplacian Graph
 Kernel for computing the gram matrix of a collection of graphs.
 
 
 Copyright (C) 2016 Imre Risi Kondor, Horace Pan
 Copyright (C) 2016 Imre Risi Kondor, Nedelina Teneva, Pramod K Mudrakarta
 
 
 The following code is a derivative work of the code from the pMMF library(https://github.com/risi-kondor/pMMF) 
 which is licensed under the GNU Public License, version 3. 
 This code therefore is also licensed under the terms of the GNU Public License, version 3. 
 ----------------------------------------------------------------------------- */

#ifndef _ThreadManager
#define _ThreadManager

#include <list>
#include <queue>
#include <mutex>

class ThreadBank;

using namespace std;


class ThreadManager{
public:

  ThreadManager(const int _maxthreads):maxthreads(_maxthreads),nthreads(0){}
  ~ThreadManager(){}

public:
  
  void enqueue(ThreadBank* bank);
  void release(ThreadBank* bank);

  int get_nthreads(){lock_guard<mutex> lock(mx); return nthreads;}

private:

  bool is_runnable(ThreadBank* bank);
  void launch(ThreadBank* bank);

public:

  int maxthreads;

private:

  mutex mx;
  int nthreads;
  list<ThreadBank*> queue;

};


#endif
