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


#ifndef _Remap
#define _Remap

#include "pMMFbase.hpp"

class Remap{
public:

  Remap(){forward=NULL; backward=NULL; n=0;}
  Remap(const Remap& x);
  Remap(Remap&& x);
  Remap& operator=(const Remap& x);
  Remap operator=(Remap&& x);
  ~Remap(){delete[] forward; delete[] backward;}

public:

  Remap(const int _n): n(_n){
    forward=new int[n]; for(int i=0; i<n; i++) forward[i]=i;
    backward=new int[n]; for(int i=0; i<n; i++) backward[i]=i;
  }

  // Remap(const int _n, const Random random);

public:

  static Remap Random(const int n);

public:

  

  int operator()(const int i) const{return forward[i];}

  Remap& swap(const int i, const int j){
    int t=forward[i]; forward[i]=forward[j]; forward[j]=t; 
    backward[forward[i]]=i; backward[forward[j]]=j;
    return *this;}

  void fixBackwardMap(){
    for(int i=0; i<n; i++) backward[forward[i]]=i;}

  string str() const;

public:
  
  int n;
  int* forward;
  int* backward;

};


#endif 
