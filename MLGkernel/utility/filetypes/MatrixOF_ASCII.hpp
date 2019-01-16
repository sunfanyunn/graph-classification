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
#ifndef _MatrixOF_ASCII
#define _MatrixOF_ASCII

#include "MatrixOF.hpp"



class MatrixOF_ASCII: public MatrixOF{
public:

  class Dense;
  class Sparse;

};



class MatrixOF_ASCII::Dense: public MatrixOF_ASCII{
public:

  Dense(const string filename, const int _nrows, const int _ncols){
    nrows=_nrows; ncols=_ncols; sparse=0; ofs.open(filename); i=0; j=0;}

public:
  
  MatrixOF& operator<<(const FIELD& v){
    ofs<<v; 
    if(++i<ncols) ofs<<" ";
    else {ofs<<"\n"; i=0; j++;}
    return *this;
  }

public:

  int i=0;
  int j=0;

};



class MatrixOF_ASCII::Sparse: public MatrixOF_ASCII{
public:

  Sparse(const string filename, const int _nrows, const int _ncols){
    nrows=_nrows; ncols=_ncols; sparse=1; ofs.open(filename);}

public:

  MatrixOF& operator<<(const IndexValueTriple& t){
    ofs<<t.i<<" "<<t.j<<" "<<t.value<<"\n";
    return *this;
  }

public:

};


#endif
