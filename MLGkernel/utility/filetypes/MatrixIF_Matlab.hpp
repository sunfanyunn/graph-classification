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
#ifndef _MatrixIF_Matlab
#define _MatrixIF_Matlab

#include "MatrixIF.hpp"
#include <matio.h>


class MatrixIF_Matlab: public MatrixIF{
public:

  class Dense;
  class Sparse;

};



class MatrixIF_Matlab::Dense: public MatrixIF_Matlab{
public:

  Dense(const string filename){
    sparse=0;
    matfile=Mat_Open(filename.c_str(),MAT_ACC_RDONLY);
    if(matfile==NULL){cout<<"Error: file cannot be opened"<<endl; return;}
    matvar=Mat_VarReadNext(matfile);
    if(matvar==NULL){cout<<"Error: no variables"<<endl; return;}
    nrows=matvar->dims[0];
    ncols=matvar->dims[1];
    next=reinterpret_cast<double*>(matvar->data);

    //need to swap to 'transpose' the array since matIO reads it in column major order
    /*int n = 5; //sqrt(nrows*ncols); 
      cout<<"size"<<n<<endl;
      for(int y = 0; y < n; ++y)
      for(int x = y+1; x < n; ++x){
      swap(next[x*n + y], next[y*n + x]); 
      }*/
      int count= nrows*ncols;
      for (int x= 0; x<nrows; ++x){
        int count_adjustment= nrows - x - 1;
        for (int y= 0, step= 1; y<ncols; ++y, step+= count_adjustment){
         int last= count - (y+x*ncols);
         int first= last - step;
         std::rotate(next + first, next + first + 1, next + last);
       }
     }

   }



 public:

  MatrixIF& operator>>(FIELD& dest){
    dest = *next;
    next++;
    return *this;
  }

public:

  mat_t* matfile;
  matvar_t* matvar;
  double *next;

};



class MatrixIF_Matlab::Sparse: public MatrixIF_Matlab{
public:

  Sparse(const string filename){
    sparse=1;
    matfile=Mat_Open(filename.c_str(),MAT_ACC_RDONLY);
    if(matfile==NULL){cout<<"Error: file "<<filename<<" cannot be opened"<<endl; return;}
    matvar=Mat_VarReadNext(matfile);
    if(matvar==NULL){cout<<"Error: no variables"<<endl; return;}
    mat_sparse_t *sparse;
    if (matvar->class_type == MAT_C_SPARSE){
      sparse = (mat_sparse_t*)matvar->data;
    }
    nrows=matvar->dims[0];
    ncols=matvar->dims[1];
    next=reinterpret_cast<double*>(sparse->data); 
    //Mat_VarPrint(matvar,1);
    //cout<<"printed"<<endl;
    Jc= sparse->jc;
    Ir = sparse->ir; 
    njc= sparse->njc; 
    ndata=sparse->ndata; 
    cout<<"read"<<endl;
  }


public:

  MatrixIF& operator>>(IndexValueTriple& dest){
    int i = indIr; int j = indJc; int c = 0;  
    for (; i < njc-1; i++ ) { 
      c= 0;
      for (; j<Jc[i+1] && j<ndata;j++) {
       c++;
       dest.i=Ir[j]; dest.j=i; 
       indJc = j+1; indIr = i; 
       break; 
     }
     if (c>0){
       break; 
     }
     indIr = i; 
   }
   dest.value =  *next++; 
   if(!dest.value || i>=njc-1){dest.i=-1; return *this;}
   return *this;
 }


public:

  mat_t* matfile;
  matvar_t* matvar;
  double *next;

  int indIr= 0; //sparse->ir; 
  int indJc =0; 
  int* Ir; 
  int* Jc; //sparse->jc;
  int njc; 
  int ndata; 

};


#endif
