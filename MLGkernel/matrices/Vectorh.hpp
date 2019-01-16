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


#ifndef _Vectorh
#define _Vectorh

#include "SparseVector.hpp"
#include <unordered_map>
#include "Cvector.hpp"


class Vectorh: public SparseVector, public unordered_map<INDEX,FIELD>, public Serializable{
public:

  Vectorh copy(){Vectorh v(n); for(auto& p:*this) v[p.first]=p.second; return v;}


public:

  Vectorh(const int _n): SparseVector(_n){}
  // Vectorh(const int _n, const Uninitialized& dummy): Vectorh(_n){};
  // Vectorh(const int _n, const Random& dummy); // obsolete 


public: // named constructors

  static Vectorh Random(const int _n, const FIELD p=0.5);
  static Vectorh Zero(const int _n){return Vectorh(_n);}
  static Vectorh Remap(const Vectorh& x, const class Remap& remap) {return Vectorh(x,remap,false);}
  static Vectorh InverseRemap(const Vectorh& x, const class Remap& remap) {return Vectorh(x,remap,true);}


public: // converters

  Vectorh(const Cvector& x);
  Vectorh(const Vectorv& x);
  Vectorh(const Vectorl& x);
  Vectorh(const Vectorl& x, const class Remap& remap, const bool inverse=false);


public: // element access

  FIELD& operator()(const int i) {assert(i<n); return (*this)[i];}

  FIELD operator()(const int i) const{assert(i<n);
    const_iterator it=find(i); if(it==end()) return 0; else return it->second;}

  FIELD read(const int i) const{assert(i<n);
    const_iterator it=find(i); if(it==end()) return 0; else return it->second;}

  bool isFilled(const int i)const {assert(i<n); const_iterator it=find(i); return (it!=end());}
  int nFilled() const {return size();}


public: // sparse methods 

  FIELD* findptr(const int i){
    iterator it=find(i); 
    if(it!=end()) return &it->second;
    return &dummyZero;
  }

  void insert(const int i, const FIELD value){assert(i<n); (*this)[i]=value;}
  void append(const int i, const FIELD value){assert(i<n); (*this)[i]=value;}
  void zero(const int i){assert(i<n); erase(i);}

  void sort(){}
  void tidy() {for(auto it=begin(); it!=end(); it++) if(it->second==0) it=erase(it);}

  void (foreach)(std::function<void(const INDEX, FIELD&)> lambda){
    for(auto& p:*this) lambda(p.first,p.second);}

  void (foreach)(std::function<void(const INDEX, const FIELD)> lambda) const{
    for(auto& p:*this) lambda(p.first,p.second);}

  void for_each(std::function<void(const INDEX, FIELD&)> lambda){
    for(auto& p:*this) lambda(p.first,p.second);}

  void for_each(std::function<void(const INDEX, const FIELD)> lambda) const{
    for(auto& p:*this) lambda(p.first,p.second);}


public: // scalar operations 

  int nnz() const {return size();}

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


  FIELD diff2(const Vectorh& x) const{
    assert(x.n==n); FIELD t=0;
    for(auto& xp:x) if(xp.second!=0) t+=(xp.second-read(xp.first))*(xp.second-read(xp.first));
    for(auto& p:*this) if(x.read(p.first)==0) t+=p.second*p.second;
    return t;
  }

  FIELD dot(const Vectorh& x) const{
    FIELD t=0; for(auto& p:x) t+=p.second*this->read(p.first); return t;}


public: // in-place operations 

  Vectorh& operator*=(const FIELD c){
    for(auto& p:*this) p.second*=c; return *this;}

  Vectorh& operator/=(const FIELD c){
    for(auto& p:*this) p.second/=c; return *this;}

  Vectorh& operator*=(const Cvector& x){
    assert(x.n==n); for(auto& p:*this) p.second*=x(p.first); return *this;}
  
  Vectorh& operator/=(const Cvector& x){
    assert(x.n==n); for(auto& p:*this) p.second/=x(p.first); return *this;}
  
  Vectorh& add(const Vectorh& x, const FIELD c=1){
    for(auto& xp:x) (*this)[xp.first]+=c*xp.second; return *this;}

  Vectorh& operator+=(const Vectorh& x){
    for(auto& xp:x) (*this)[xp.first]+=xp.second; return *this;}
  Vectorh& operator-=(const Vectorh& x){
    for(auto& xp:x) (*this)[xp.first]+=xp.second; return *this;}
  Vectorh& operator*=(const Vectorh& x){
    for(auto& xp:x) (*this)[xp.first]*=xp.second; return *this;}
  Vectorh& operator/=(const Vectorh& x){
    for(auto& xp:x) (*this)[xp.first]/=xp.second; return *this;}

public:

  string str(const Sparse dummy) const;
  string str(const Dense dummy) const {return Vector::str(Dense());}
  string str() const{return str(Sparse());}

  static string classname();
  Vectorh(Bifstream& ifs);
  void serialize(Bofstream& ofs) const;
  void serialize(Rstream& rstream) const;

public:  

  static FIELD dummyZero;

};

#endif
