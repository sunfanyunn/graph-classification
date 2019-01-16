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


#include "Vector.hpp"
#include <vector>


string Vector::str(const Dense dummy) const{
  ostringstream stream;
  for(int i=0; i<n; i++) stream<<(*this)(i)<<endl; 
  return stream.str();  
}


string Vector::str(const Sparse dummy) const{
  ostringstream stream;
  for(int i=0; i<n; i++) if((*this)(i)!=0) stream<<(*this)(i)<<endl; 
  return stream.str();  
}


ostream& operator<<(ostream& stream, const Vector& x) {stream<<x.str(); return stream;}
//  {stream<<"("; for (int i=0; i<x.n-1; i++) stream<<x(i)<<","; stream<<x(x.n-1)<<")"; return stream;}

