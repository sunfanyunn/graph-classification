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


#ifndef _Vector
#define _Vector

#include "pMMFbase.hpp"
#include <functional>

class Vector{ //: public Serializable{
public:

  Vector(const int _n): n(_n){}

public:

  virtual FIELD& operator()(const int n)=0;
  virtual FIELD operator()(const int n) const=0;
  virtual FIELD read(const int i) const {return (*this)(i);}

  //virtual void (foreach)(std::function<void(const INDEX, FIELD&)> lambda)=0;
  //virtual void (foreach)(std::function<void(const INDEX, const FIELD)> lambda) const=0;

  virtual bool isFilled (const int i)const =0;
  virtual int nFilled() const=0;

public:

  virtual int nnz() const=0;

  virtual int argmax() const=0;
  virtual int argmax_abs() const=0;

  virtual FIELD norm2() const=0;
  // FIELD diff2(const VECTOR& x)=0;

public:

  //virtual void serialize(Bofstream& ofs) const=0;
  //virtual void serialize(Rstream& rstream) const=0;

  virtual string str(const Dense dummy) const; 
  virtual string str(const Sparse dummy) const; 
  virtual string str() const{return str(Dense());};

public:  

  int n;

};


ostream& operator<<(ostream& stream, const Vector& x);




#endif
