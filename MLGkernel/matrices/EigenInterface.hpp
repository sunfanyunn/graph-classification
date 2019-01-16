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


#ifndef _EigenInterface
#define _EigenInterface

#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <unsupported/Eigen/MatrixFunctions>
	
// The purpose of these adaptors is to avoid having to include Eigen/Dense or Eigen/Core in any of the 
// header files of the native vector/matrix classes, which would slow down compilation.  

typedef Eigen::SparseMatrix<FIELD> EigenSparseMatrix;

class EigenVectorXdAdaptor: public Eigen::VectorXd{
public:
  EigenVectorXdAdaptor(const Eigen::VectorXd& M): Eigen::VectorXd(M){}
};

class EigenMatrixXdAdaptor: public Eigen::MatrixXd{
public:
  EigenMatrixXdAdaptor(const Eigen::MatrixXd& M): Eigen::MatrixXd(M){}
};


#endif
