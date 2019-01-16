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


#ifndef _AtomicCmatrix
#define _AtomicCmatrix

#include "Cmatrix.hpp"

class AtomicCmatrix: public Cmatrix{
public:

  using Cmatrix::Cmatrix;

  //Cmatrix(const Cmatrix& x); // these will need to be implemented 
  //Cmatrix(Cmatrix&& x);
  //Cmatrix& operator=(const Cmatrix& x);
  //Cmatrix& operator=(Cmatrix&& x);
  //~Cmatrix(){delete[] array;}
  
public:

public: // in-place operations 

  Cmatrix& operator+=(const Cmatrix& x){
    lock_guard<mutex> lock(mx);
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    for(int i=0; i<nrows*ncols; i++) array[i]+=x.array[i]; return *this;}


public:

  mutex mx;

};


#endif
