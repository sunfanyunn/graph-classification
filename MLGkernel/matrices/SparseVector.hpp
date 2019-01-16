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


#ifndef _SparseVector
#define _SparseVector

#include "Vector.hpp"



class SparseVector: public Vector{
public:

  SparseVector(const int _n): Vector(_n){}


public:

  virtual FIELD& operator()(const int i)=0;
  virtual FIELD operator()(const int i) const=0;
  virtual FIELD read(const int i) const=0;

  virtual void insert(const int i, const FIELD value)=0;
  virtual void append(const int i, const FIELD value)=0;
  virtual void zero(const int i)=0;

  virtual void sort()=0;
  virtual void tidy()=0;
  //virtual int nnz() const=0;
};






#endif
