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


#include "MatrixX.hpp"


/*
template<>
void MatrixX<Vectorv>::serialize(Bofstream& ofs) const{
  ofs.tag("MatrixXv",0);
  ofs.write(nrows);
  ofs.write(ncols);
  for(int j=0; j<ncols; j++) 
    column[j]->serialize(ofs);
};

template<>
void MatrixX<Vectorl>::serialize(Bofstream& ofs) const{
  ofs.tag("MatrixXl",0);
  ofs.write(nrows);
  ofs.write(ncols);
  for(int j=0; j<ncols; j++) 
    column[j]->serialize(ofs);
};

template<>
void MatrixX<Vectorh>::serialize(Bofstream& ofs) const{
  ofs.tag("MatrixXh",0);
  ofs.write(nrows);
  ofs.write(ncols);
  for(int j=0; j<ncols; j++) 
    column[j]->serialize(ofs);
};
*/

/*
template<>
MatrixX<Vectorv>::MatrixX(Bifstream& ifs):SparseMatrix(0,0){
  ifs.check("MatrixXv",0);
  ifs.read(nrows);
  ifs.read(ncols);
  for(int j=0; j<ncols; j++)
    column.push_back(new Vectorv(ifs));
};

template<>
MatrixX<Vectorl>::MatrixX(Bifstream& ifs):SparseMatrix(0,0){
  ifs.check("MatrixXl",0);
  ifs.read(nrows);
  ifs.read(ncols);
  for(int j=0; j<ncols; j++)
    column.push_back(new Vectorl(ifs));
};

template<>
MatrixX<Vectorh>::MatrixX(Bifstream& ifs):SparseMatrix(0,0){
  ifs.check("MatrixXh",0);
  ifs.read(nrows);
  ifs.read(ncols);
  for(int j=0; j<ncols; j++)
    column.push_back(new Vectorh(ifs));
};
*/






