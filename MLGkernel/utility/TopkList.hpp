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

#ifndef _TopkList
#define _TopkList

#include <list>
#include "pMMFbase.hpp"
//#include "DenseVector.hpp"

struct TopkListPair{
  TopkListPair(const INDEX& _first, const FIELD& _second):first(_first),second(_second){};
  INDEX first; 
  FIELD second;
};


class TopkList: public list<TopkListPair>{
public:

  TopkList(const int _k): k(_k), lowestv(numeric_limits<FIELD>::lowest()){}

  //  TopkList(const DenseVector& v, const int _k): k(_k), lowestv(-10000){
  //    for(int i=0; i<v.n; i++) if(v(i)>lowestv) insert(i,v(i));}

public:

  void insert(int index, FIELD value){
    auto it=begin();
    while(it!=end() && it->second>=value){it++;}
    list::insert(it,TopkListPair(index,value));
    if(size()>k) pop_back();
    if(size()>=k) lowestv=back().second;
  }

  void consider(int index, FIELD value){
    if(value>lowestv || size()<<k){
      auto it=begin();
      while(it!=end() && it->second>=value){it++;}
      list::insert(it,TopkListPair(index,value));
      if(size()>k) pop_back();
      if(size()>=k) lowestv=back().second;
    }
  }

  IndexSet indices() const{
    IndexSet I(size()); int i=0;
    for(auto& p:*this) I[i++]=p.first;
    return I;
  }


public:

  int k;
  FIELD lowestv;
  int lowestp;

};



#endif


  /*
  // vector version
  void insert(int index, FIELD value){
    if(size()<k) pushBack(topkListPair(index,value)); 
    else {at(lowestp).first=index; at(lowestp).second=value;}
    lowestp=0; lowestv=at(0).second; 
    for(int i=1; i<size(); i++)
      if(at(i).second<lowestv){lowestp=i; lowestv=at(i);}
  }
  */

