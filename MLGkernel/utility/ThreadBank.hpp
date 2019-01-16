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

#ifndef _ThreadBank
#define _ThreadBank

#include <thread>
#include "pMMFbase.hpp"
#include "ThreadManager.hpp"

using namespace std;


//extern mutex cout_mutex;
extern ThreadManager threadManager;


class ThreadBank{
public:

  ThreadBank()=delete;
  
  ThreadBank(const int _maxthreads=1000, const int _maxprivileged=1): 
    maxthreads(_maxthreads), maxprivileged(_maxprivileged), nthreads(0), nprivileged(0) {gate.lock();}; 

  ~ThreadBank(){for(auto& th:threads) th.join();}


public:

  template<class FUNCTION, class OBJ>
  void add(FUNCTION lambda, const OBJ x){
    lock_guard<mutex> lock(mx); //                                   unnecessary if called from a single thread
    threadManager.enqueue(this);
    gate.lock(); //                                                  gate can only be unlocked by threadManager
    nthreads++;
    threads.push_back(thread([this,lambda](OBJ _x){lambda(_x); nthreads--; threadManager.release(this);},x));
    #ifdef _THREADBANKVERBOSE
    printinfo();
    #endif
  }


  template<class FUNCTION, class OBJ1, class OBJ2>
  void add(FUNCTION lambda, const OBJ1 x1, const OBJ2 x2){
    lock_guard<mutex> lock(mx);
    threadManager.enqueue(this);
    gate.lock();
    nthreads++;
    threads.push_back(thread([this,lambda](OBJ1 _x1, OBJ2 _x2){
			       lambda(_x1,_x2); nthreads--; threadManager.release(this);},x1,x2));
    #ifdef _THREADBANKVERBOSE
    printinfo();
    #endif
  }


  template<class FUNCTION, class OBJ1, class OBJ2, class OBJ3>
  void add(FUNCTION lambda, const OBJ1 x1, const OBJ2 x2, const OBJ3 x3){
    lock_guard<mutex> lock(mx);
    threadManager.enqueue(this);
    gate.lock();
    nthreads++;
    threads.push_back(thread([this,lambda](OBJ1 _x1, OBJ2 _x2, OBJ3 _x3){
			       lambda(_x1,_x2,_x3); nthreads--; threadManager.release(this);},x1,x2,x3));
    #ifdef _THREADBANKVERBOSE
    printinfo();
    #endif
  }


  bool is_ready(){return nthreads<maxthreads;}


  void printinfo(){
    CoutLock lock;
    cout<<"    (threads: "<<nthreads-nprivileged<<"+"<<nprivileged<<" local, ";
    cout<<threadManager.get_nthreads()<<" global)"<<endl;
  }

  
public:

  mutex mx;
  mutex gate;
  atomic<int> nthreads;
  int nprivileged=0; //                                               only to be touched by threadManager
  int maxthreads=4;
  int maxprivileged=1;

  vector<thread> threads;

};





#endif
