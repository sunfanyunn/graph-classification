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
#ifndef _MatrixIF
#define _MatrixIF

#include "pMMFbase.hpp"
#include <fstream>

class MatrixIF{
public:
  
  ~MatrixIF(){ifs.close();}

public:

  virtual void rewind(){}

  virtual MatrixIF& operator>>(FIELD& v){
    cout<<"Error: operator>>(FIELD& ) not supported in sparse matrix input files."<<endl; 
    return *this;}

  virtual MatrixIF& operator>>(IndexValueTriple& t){
    cout<<"Error: operator>>(IndexValueTriple& ) not supported in dense matrix input files."<<endl; 
    return *this;}


public:

  bool sparse;
  ifstream ifs;
  int nrows;
  int ncols;

};


#endif 
