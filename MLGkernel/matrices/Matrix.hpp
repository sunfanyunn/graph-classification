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


#ifndef _Matrix
#define _Matrix

#include "pMMFbase.hpp"
#include <functional>

class Matrix{
public:

  virtual ~Matrix(){}

public: // constructors

  Matrix(const int _nrows, const int _ncols): nrows(_nrows), ncols(_ncols) {}


public: // member access
  
  virtual FIELD& operator()(const int i, const int j)=0; 
  virtual FIELD operator()(const int i, const int j) const=0; 
  virtual FIELD read(const int i, const int j) const=0; 
  virtual bool isFilled(const int i, const int j) const=0; 
  virtual int nFilled() const=0;
  virtual bool isSparse() const=0;

  //virtual void (foreach)(std::function<void(const INDEX, const INDEX,FIELD&)> lambda)=0;
  //virtual void (foreach)(std::function<void(const INDEX, const INDEX, const FIELD)> lambda) const=0;
  virtual void foreach_in_column(const int j, std::function<void(const INDEX, FIELD&)> lambda)=0;
  virtual void foreach_in_column(const int j, std::function<void(const INDEX, const FIELD)> lambda) const=0;

public: // scalar valued operations 

  virtual int nnz() const=0;

public:

  virtual Cmatrix dot(const Matrix& x) const; // {};

public:

  virtual void saveTo(MatrixOF& file) const=0;

  virtual string str(const Dense dummy) const;
  virtual string str(const Sparse dummy) const;
  virtual string str() const{return str(Dense());}

   
  
public:

  int nrows;
  int ncols;

};


ostream& operator<<(ostream& stream, const Matrix& x); 



class SparseMatrix: public Matrix{
public:
  using Matrix::Matrix;
  bool isSparse() const {return true;}
};



class DenseMatrix: public Matrix{
public:
  using Matrix::Matrix;
  bool isFilled(const int i, const int j) const {return true;}
  int nFilled() const {return nrows*ncols;}
  bool isSparse() const {return false;}
};




// virtual Matrix* newof()=0;

//virtual Cmatrix Cmatrix() const=0;
// virtual SparseMatrixX<SparseVectorv> MatrixXv() const=0;
// virtual SparseMatrixX<SparseVectorl> MatrixXl() const=0;
// virtual SparseMatrixX<SparseVectorh> MatrixXh() const=0;



#endif 
