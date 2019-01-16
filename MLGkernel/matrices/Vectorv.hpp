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


#ifndef _Vectorv
#define _Vectorv

#include "SparseVector.hpp"
#include "Cvector.hpp"
#include <algorithm>

class Vectorv: public SparseVector, public vector<SVpair>, Serializable{
public:

  //Vectorv(Vectorv&& x): SparseVector(x.n), sorted(x.sorted){vector<SVpair>=std::move(x);}
  Vectorv copy(){Vectorv v(n); v.resize(size()); for(int i=0; i<size(); i++) v[i]=(*this)[i]; return v;}


public: // constructors

  Vectorv(const int _n): SparseVector(_n), sorted(1){}
  // Vectorv(const int _n, const Random& dummy); // obsolete


public: // named constructors

  static Vectorv Zero(const int _n){return Vectorv(_n);}
  static Vectorv Filled(const int _n, const FIELD t){
    Vectorv v(_n); v.resize(_n); for(int i=0; i<_n; i++) v[i]=SVpair(i,t); return v;}
  static Vectorv Random(const int _n, const FIELD p=0.5);
  static Vectorv Remap(const Vectorv& x, const class Remap& remap) {return Vectorv(x,remap,false);}
  static Vectorv InverseRemap(const Vectorv& x, const class Remap& remap) {return Vectorv(x,remap,true);}


public: // converters 

  Vectorv(const Cvector& x);
  Vectorv(const Vectorl& x);
  Vectorv(const Vectorh& x);
  Vectorv(const Vectorv& x, const class Remap& remap, const bool inverse=false);


public: // element access 

  FIELD& operator()(const int i) {assert(i<n);
    auto it=begin(); while(it!=end() && it->first!=i){it++;}
    if(it==end()) {push_back(SVpair(i,0)); sorted=0; return back().second;}
    return it->second;
  }

  FIELD operator()(const int i) const {assert(i<n);
    auto it=begin(); while(it!=end() && it->first!=i){it++;}
    if(it==end()) return 0; else return it->second;
  }

  FIELD read(const int i) const {assert(i<n);
    auto it=begin(); while(it!=end() && it->first!=i){it++;}
    if(it==end()) return 0; else return it->second;
  }

  void (foreach)(std::function<void(const INDEX, FIELD&)> lambda){
    for(auto& p:*this) lambda(p.first,p.second);}

  void (foreach)(std::function<void(const INDEX, const FIELD)> lambda) const{
    for(auto& p:*this) lambda(p.first,p.second);}

  void for_each(std::function<void(const INDEX, FIELD&)> lambda){
    for(auto& p:*this) lambda(p.first,p.second);}

  void for_each(std::function<void(const INDEX, const FIELD)> lambda) const{
    for(auto& p:*this) lambda(p.first,p.second);}

  bool isFilled(const int i) const {assert(i<n);
    auto it=begin(); while(it!=end() && it->first!=i){it++;}
    return (it!=end());}

  int nFilled() const {return size();}


public: // sparse vector methods

  FIELD* findptr(const int i){
    auto it=begin(); 
    if(sorted) {while(it!=end() && it->first<i) {it++;} if(it!=end() && it->first==i) return &it->second;}
    else while(it!=end()) {if(it->first==i) return &it->second; it++;}
    return &dummyZero;
  }

  void insert(const int i, const FIELD value) {assert(i<n); push_back(SVpair(i,value)); sorted=0;}

  void append(const int i, const FIELD value) {assert(i<n); push_back(SVpair(i,value));}

  void zero(const int i){assert(i<n);
    auto it=begin(); while(it!=end() && it->first!=i){it++;}
    if(it!=end()) it->second=0;}

  void sort() {if(!sorted) std::sort(begin(), end(), SVpairComparator); sorted=1;}

  void tidy(){
    sort(); auto v(*this); clear(); int i=-1;
    for(auto& p:v) if(p.first!=i && p.second!=0) {push_back(p); i=p.first;}}
  

public: // scalar-valued operations 

  int nnz() const {const_cast<Vectorv&>(*this).tidy(); return size();}

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


  FIELD diff2(const Vectorv& x) const {
    assert(x.n==n); FIELD t=0; 
    const_cast<Vectorv&>(*this).sort(); 
    const_cast<Vectorv&>(x).sort(); 
    auto it=begin();
    for(auto& xv:x){
      while(it!=end() && it->first<xv.first) {t+=(it->second)*(it->second); it++;}
      if(it!=end() && it->first==xv.first) {t+=(it->second-xv.second)*(it->second-xv.second);it++;}
      else {t+=xv.second*xv.second;}
    }
    return t; 
  }

  FIELD dot(const Vectorv& x) const{
    FIELD t=0; 
    const_cast<Vectorv&>(*this).sort(); 
    const_cast<Vectorv&>(x).sort(); 
    auto xit=x.begin();
    for(auto& v:*this){
      int i=v.first;
      while(xit!=x.end() && xit->first<i) xit++;
      if(xit==x.end()) break;
      if(xit->first==i) t+=v.second*xit->second;
    }
    return t;
  }


public: // in-place operations 

  Vectorv& operator*=(const FIELD c){for(auto& p:*this) p.second*=c; return *this;}
  Vectorv& operator/=(const FIELD c){for(auto& p:*this) p.second/=c; return *this;}

  Vectorv& operator*=(const Cvector& x){assert(x.n==n); for(auto& p:*this) p.second*=x(p.first); return *this;}
  Vectorv& operator/=(const Cvector& x){assert(x.n==n); for(auto& p:*this) p.second/=x(p.first); return *this;}

  Vectorv& add(const Vectorv& x, const FIELD c=1){
    sort();
    const_cast<Vectorv&>(x).sort(); 
    vector<SVpair> newbies;
    auto it=begin();
    for(auto& xv:x){
      while(it!=end() && it->first<xv.first) it++;
      if(it!=end() && it->first==xv.first) it->second+=c*xv.second;
      else newbies.push_back(SVpair(xv.first,c*xv.second));
    }
    if(newbies.size()!=0) {vector<SVpair>::insert(end(),newbies.begin(),newbies.end()); sorted=0;}
    return *this; 
  }

  Vectorv& operator+=(const Vectorv& x){return add(x,1);}
  Vectorv& operator-=(const Vectorv& x){return add(x,-1);}

  Vectorv& operator*=(const Vectorv& x){
    sort();
    const_cast<Vectorv&>(x).sort(); 
    auto it=x.begin();
    for(auto& p:*this){
      while(it!=x.end() && it->first<p.first) it++;
      if(it!=x.end() && it->first==p.first) p.second*=it->second;
      else p.second=0;
    }
    return *this; 
  }

  Vectorv& operator/=(const Vectorv& x){
    sort();
    const_cast<Vectorv&>(x).sort(); 
    auto it=x.begin();
    for(auto& p:*this){
      while(it!=x.end() && it->first<p.first) it++;
      if(it!=x.end() && it->first==p.first) p.second/=it->second;
      else p.second=1.0/0.0;
    }
    return *this; 
    }


public: // I/O 

  string str(const Sparse dummy) const;
  string str(const Dense dummy) const {return Vector::str(Dense());}
  string str() const{return str(Sparse());}
  static string classname();

  Vectorv(Bifstream& ifs);
  void serialize(Bofstream& ofs) const;
  void serialize(Rstream& rstream) const;

public:  

  bool sorted;

  static FIELD dummyZero;

};






#endif
