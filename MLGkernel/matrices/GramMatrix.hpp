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


#ifndef _GramMatrix
#define _GramMatrix

#include <unordered_map>
#include "MatrixX.hpp"

template<class MATRIX>
class GramMatrix: public MATRIX{
public:

  using class MATRIX::MATRIX;

public:

  template<class COLTYPE>
  GramMatrix(MatrixX<COLTYPE>& A)


public:


};


template<class MATRIX>
template<class COLTYPE>
GramMatrix<MATRIX>::GramMatrix(MatrixX<COLTYPE>& A): MATRIX(MATRIX::Zero(A.nrows,A.nrows)){
  assert(A.nrows==A.ncols); // assumption: A is symmetric
  for(int i=0; i<nrows; i++){
    if(A.column[i]->nFilled>0.2*nrows){
      for(int j=0; j<=i; j++){
	(*this)(i,j)=A.column[i]->dot(*A.column[j]);
	(*this)(j,i)=(*this)(i,j);
      }
    }else{
      unordered_map neighbors;
      A.column[i]->for_each([&A,&neighbors](int j, FIELD dummy){
			      A.column[j]->for_each([&neighbors](int k, FIELD dummy){neighbors.insert(k);});
			    });
      for(auto j:neighbors){
	(*this)(i,j)=A.column[i]->dot(*A.column[j]);
	(*this)(j,i)=(*this)(i,j);
      }
    }
  }

}



#endif
