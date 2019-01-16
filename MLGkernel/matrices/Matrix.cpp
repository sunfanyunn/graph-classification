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


#include "Matrix.hpp"
#include "Cmatrix.hpp"


Cmatrix Matrix::dot(const Matrix& x) const{
  return Cmatrix(0,0);
}


string Matrix::str(const Dense dummy) const{
  ostringstream stream;
  stream.precision(3); 
  stream.setf(ios_base::fixed, ios_base::floatfield);
  for(int i=0; i<nrows; i++){stream<<"[ ";
    for(int j=0; j<ncols; j++) {stream.width(6); stream<<this->read(i,j)<<" ";}
    stream<<" ]\n";}
  return stream.str();
}


string Matrix::str(const Sparse dummy) const{
  ostringstream stream;
  for(int i=0; i<nrows; i++)
    for(int j=0; j<ncols; j++)
      if((*this)(i,j)!=0) stream<<"("<<i<<","<<j<<") : "<<(*this)(i,j)<<endl;
  return stream.str();
}


ostream& operator<<(ostream& stream, const Matrix& x){stream<<x.str(); return stream;}
