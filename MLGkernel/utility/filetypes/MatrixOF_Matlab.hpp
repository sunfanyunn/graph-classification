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
#ifndef _MatrixOF_Matlab
#define _MatrixOF_Matlab

#include "MatrixOF.hpp"
#include "matio.h"


class MatrixOF_Matlab: public MatrixOF{
public:

  using MatrixOF::MatrixOF;

  class Dense;
  class Sparse;

};



class MatrixOF_Matlab::Dense: public MatrixOF_Matlab{
public:

  Dense(const string filename, const int _nrows, const int _ncols, const FIELD* array){
    nrows=_nrows; ncols=_ncols; sparse=0; 
    matfile = Mat_CreateVer(filename.c_str(),NULL,MAT_FT_DEFAULT);
    if(matfile==NULL) { cout<<"Error: Cannot open file for write!"<<endl; return; }
    size_t dims[2]; dims[0]=nrows; dims[1]=ncols;
    matvar = Mat_VarCreate("M_dense",MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims,(void*)array,0);
    if(matvar==NULL) { cout<<"Error creating matrix variable!"<<endl; return; }
    Mat_VarWrite(matfile,matvar,MAT_COMPRESSION_NONE);
    Mat_VarFree(matvar);
    Mat_Close(matfile);
  }

  ~Dense(){}

public:

  MatrixOF& operator<<(const FIELD& v){
    return *this;
  }

public:
  mat_t *matfile;
  matvar_t* matvar;

};



class MatrixOF_Matlab::Sparse: public MatrixOF_Matlab{
public:

  Sparse(const string filename, const CSCmatrix &cscmatrix){
    nrows=cscmatrix.nrows; ncols=cscmatrix.ncols; sparse=1; 
    size_t dims[2]; dims[0]=nrows; dims[1]=ncols;
    mat_sparse_t sparsemat = {0,};
    sparsemat.nzmax = cscmatrix.nnz;
    sparsemat.nir = cscmatrix.nnz;
    sparsemat.ir = cscmatrix.ir;
    sparsemat.njc = cscmatrix.ncols+1;
    sparsemat.jc = cscmatrix.jc;
    sparsemat.ndata = cscmatrix.nnz;
    sparsemat.data = cscmatrix.val;

    matfile = Mat_CreateVer(filename.c_str(),NULL,MAT_FT_DEFAULT);
    if(matfile==NULL) { cout<<"Error: Cannot open file for write!"<<endl; return; }
    matvar = Mat_VarCreate("M_sparse",MAT_C_SPARSE,MAT_T_DOUBLE,2,dims,&sparsemat,0);
    if(matvar==NULL) { cout<<"Error creating matrix variable!"<<endl; return; }
    Mat_VarWrite(matfile,matvar,MAT_COMPRESSION_NONE);
    Mat_VarFree(matvar);
    Mat_Close(matfile);
  }

  ~Sparse() {}

public:

  MatrixOF& operator<<(const IndexValueTriple& t){
    return *this;
  }

public:
  mat_t* matfile;
  matvar_t* matvar;


};



#endif
