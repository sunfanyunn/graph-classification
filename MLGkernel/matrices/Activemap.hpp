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


#ifndef _Activemap
#define _Activemap

#include "Remap.hpp"

class Activemap: public Remap{
public:

  Activemap(const int n=1): Remap(n), nactive(n){}

public:

  int random();

  void remove(const int i){
    if(backward[i]!=nactive-1) swap(backward[i],nactive-1); 
    nactive--; 
   }

  bool isactive(const int i) const {return(backward[i]<nactive);}

public:

  int nactive=0;
  
};


#endif
