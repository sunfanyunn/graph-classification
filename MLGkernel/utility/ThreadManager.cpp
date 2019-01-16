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

#include "ThreadManager.hpp"
#include "ThreadBank.hpp"


void ThreadManager::enqueue(ThreadBank* bank){
  lock_guard<mutex> lock(mx);
  if(is_runnable(bank)) launch(bank);
  else queue.push_back(bank);
}


void ThreadManager::release(ThreadBank* bank){
  lock_guard<mutex> lock(mx);
  if(bank->nprivileged>0) bank->nprivileged--;
  else nthreads--;
  for(auto it=queue.begin(); it!=queue.end(); it++)
    if(is_runnable(*it)){
      launch(*it);
      it=queue.erase(it);
    }
  //  auto it=find_if(queue.begin(),queue.end(),[this](ThreadBank* bank){return is_runnable(bank);});
  // if(it==queue.end()) return;
  // ThreadBank* bank=*it;
  // queue.erase(it);
  // launch(bank);
}


bool ThreadManager::is_runnable(ThreadBank* bank){
  return bank->is_ready() && (bank->nprivileged<bank->maxprivileged || nthreads<maxthreads) ;
}


void ThreadManager::launch(ThreadBank* bank){
  if(bank->nprivileged<bank->maxprivileged) bank->nprivileged++;
  else nthreads++;
  bank->gate.unlock();
}


  /*
  void addBank(const ThreadBank* bank){
    lock_guard<mutex> lock(mx);
    banks.push_front(bank);
  }

  void removeBank(const ThreadBank* bank){
    lock_guard<mutex> lock(mx);
    banks.remove(bank);
  }
  */
