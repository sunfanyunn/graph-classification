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


#ifndef _Vectorl
#define _Vectorl

#include "SparseVector.hpp"
#include <list>
#include "Cvector.hpp"


class Vectorl: public SparseVector, public list<SVpair>, public Serializable{
public:

  Vectorl copy(){Vectorl v(n); v.resize(size()); for(auto& p:*this) v.push_back(p); return v;}

public:

  Vectorl(const int _n): SparseVector(_n){}
  // Vectorl(const int _n, const Uninitialized& dummy): Vectorl(_n){};
  // Vectorl(const int _n, const Random& dummy); // deprecated 


public: // named constructors

  static Vectorl Zero(const int _n){return Vectorl(_n);}
  static Vectorl Filled(const int _n, const FIELD t){
    Vectorl v(_n); for(int i=0; i<_n; i++) v.push_back(SVpair(i,t)); return v;}
  static Vectorl Random(const int _n, const FIELD p=0.5);
  static Vectorl Remap(const Vectorl& x, const class Remap& remap) {return Vectorl(x,remap,false);}
  static Vectorl InverseRemap(const Vectorl& x, const class Remap& remap) {return Vectorl(x,remap,true);}



public: // converters

  Vectorl(const Cvector& x);
  Vectorl(const Vectorv& x);
  Vectorl(const Vectorh& x);
  Vectorl(const Vectorl& x, const class Remap& remap, const bool inverse=false);


public: // element access 

  FIELD& operator()(const int i){assert(i<n);
    auto it=begin(); while(it!=end() && it->first<i){it++;} // list is always ordered 
    if(it==end() || it->first>i) it=list<SVpair>::insert(it,SVpair(i,0));
    return it->second;
  }

  FIELD operator()(const int i) const{assert(i<n);
    auto it=begin(); while(it!=end() && it->first<i){it++;}
    if(it==end() || it->first>i) return 0; else return it->second;
  }

  FIELD read(const int i) const{assert(i<n);
    auto it=begin(); while(it!=end() && it->first<i){it++;}
    if(it==end() || it->first>i) return 0; else return it->second;
  }

  bool isFilled(const int i)const {assert(i<n);
    auto it=begin(); while(it!=end() && it->first<i){it++;}
    return (it!=end() && it->first==i);}

  int nFilled() const {return size();}

  void (foreach)(std::function<void(const INDEX, FIELD&)> lambda){
    for(auto& p:*this) lambda(p.first,p.second);}

  void (foreach)(std::function<void(const INDEX, const FIELD)> lambda) const{
    for(auto& p:*this) lambda(p.first,p.second);}

  void for_each(std::function<void(const INDEX, FIELD&)> lambda){
    for(auto& p:*this) lambda(p.first,p.second);}

  void for_each(std::function<void(const INDEX, const FIELD)> lambda) const{
    for(auto& p:*this) lambda(p.first,p.second);}


public: // sparse methods

  FIELD* findptr(const int i){
    auto it=begin(); while(it!=end() && it->first<i) {it++;} 
    if(it!=end() && it->first==i) return &it->second;
    return &dummyZero;
  }

  void insert(const int i, const FIELD value) {assert(i<n); (*this)(i)=value;}

  void append(const int i, const FIELD value) {assert(i<n); push_back(SVpair(i,value));}

  void zero(const int i){assert(i<n); 
    auto it=begin(); while(it!=end() && it->first<i){it++;}
    if(it!=end() && it->first==i) erase(it);}

  void sort(){}

  void tidy() {int i=-1;
    for(auto it=begin(); it!=end(); it++){
      while(it->first==i || it->second==0) it=erase(it); 
      if(it!=end()) i=it->first;}
  }


public: // scalar-valued operations 

  int nnz() const {const_cast<Vectorl&>(*this).tidy(); return size();}

  int argmax() const{
    if(size()==0) return 0;
    int best=begin()->first; FIELD max=begin()->second;
    for(auto p:*this) if(p.second>max) {best=p.first; max=p.second;}
    return best;
  }

  int argmax_abs() const{
    if(size()==0) return 0;
    int best=begin()->first; FIELD max=fabs(begin()->second);
    for(auto p:*this) if(fabs(p.second)>max) {best=p.first; max=fabs(p.second);}
    return best;
  }

  FIELD norm2() const{
    FIELD t=0; for(auto& p:*this) t+=p.second*p.second; return t;}

  FIELD norm() const { return sqrt(norm2()); }


  FIELD diff2(const Vectorl& x) const {
    assert(x.n==n); FIELD t=0;
    auto it=begin();
    for(auto& xv:x){
      while(it!=end() && it->first<xv.first) {t+=(it->second)*(it->second); it++;}
      if(it!=end() && it->first==xv.first) {t+=(it->second-xv.second)*(it->second-xv.second); it++;}
      else t+=xv.second*xv.second;
    }
    return t; 
  }

  FIELD dot(const Vectorl& x) const{
    FIELD t=0; 
    auto xit=x.begin();
    for(auto& v:*this){
      int i=v.first;
      while(xit!=x.end() && xit->first<i) xit++;
      if(xit==x.end()) break;
      if(xit->first==i) t+=v.second*xit->second;
    }
    return t;
  }

  FIELD dot(const Cvector& x) const{
    FIELD t=0; 
    for(auto& v:*this){
      t+=v.second*x(v.first);;
    }
    return t;
  }


public: // in-place operations 

  Vectorl& operator*=(const FIELD c){
    for(auto& p:*this) p.second*=c; return *this;}

  Vectorl& operator/=(const FIELD c){
    for(auto& p:*this) p.second/=c; return *this;}

  Vectorl& operator*=(const Cvector& x){
    assert(x.n==n); for(auto& p:*this) p.second*=x(p.first); return *this;}

  Vectorl& operator/=(const Cvector& x){
    assert(x.n==n); for(auto& p:*this) p.second/=x(p.first); return *this;}

  Vectorl& add(const Vectorl& x, const FIELD c=1){
    auto it=begin();
    for(auto& xv:x){
      while(it!=end() && it->first<xv.first) it++;
      if(it!=end() && it->first==xv.first) it->second+=c*xv.second;
      else {list<SVpair>::insert(it,SVpair(xv.first,c*xv.second));}
    }
    return *this;
  }

  Vectorl& operator+=(const Vectorl& x){return add(x,1);}
  Vectorl& operator-=(const Vectorl& x){return add(x,-1);}

  Vectorl& operator*=(const Vectorl& x){
    auto it=x.begin();
    for(auto& p:*this){
      while(it!=x.end() && it->first<p.first) it++;
      if(it!=x.end() && it->first==p.first) p.second*=it->second;
      else p.second=0;
    }
    return *this; 
  }

  Vectorl& operator/=(const Vectorl& x){
    auto it=x.begin();
    for(auto& p:*this){
      while(it!=x.end() && it->first<p.first) it++;
      if(it!=x.end() && it->first==p.first) p.second/=it->second;
      else p.second=1.0/0.0;
    }
    return *this; 
    }


public:

  string str(const Sparse dummy) const;
  string str(const Dense dummy) const {return Vector::str(Dense());}
  string str() const{return str(Sparse());}

  static string classname();
  Vectorl(Bifstream& ifs);
  void serialize(Bofstream& ofs) const;
  void serialize(Rstream& rstream) const;

public:  

  static FIELD dummyZero; 

};


#endif
