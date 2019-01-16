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
#ifndef _MatrixOF
#define _MatrixOF

#include "pMMFbase.hpp"
#include <fstream>

class MatrixOF{
public:
  
  //MatrixOF(const char* filename, const int _nrows, const int _ncols): 
  //  nrows(_nrows), ncols(_ncols), ofs(filename){}

  ~MatrixOF(){ofs.close();}

public:

  virtual MatrixOF& operator<<(const FIELD& v){
    cout<<"Error: operator<<(FIELD& ) not supported in sparse matrix output files."<<endl; 
    return *this;
  };

  virtual MatrixOF& operator<<(const IndexValueTriple& t){
    cout<<"Error: operator<<(IndexValueTriple& ) not supported in dense matrix output files."<<endl; 
    return *this;
  };

public:

  bool sparse=0;
  ofstream ofs;
  int nrows;
  int ncols;

};


#endif 
