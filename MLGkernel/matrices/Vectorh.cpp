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


#include "Vectorv.hpp"
#include "Vectorl.hpp"
#include "Vectorh.hpp"
#include "Cvector.hpp"
#include <random>

extern default_random_engine randomNumberGenerator;


/*
Vectorh::Vectorh(const int _n, const class Random& dummy): SparseVector(_n){
  uniform_real_distribution<FIELD> distr(0,1);
  for(int i=0; i<n; i++) 
    if(distr(randomNumberGenerator)<=dummy.p) (*this)[i]=distr(randomNumberGenerator);
}
*/

Vectorh Vectorh::Random(const int _n, const FIELD p){
  Vectorh v(_n); 
  uniform_real_distribution<FIELD> distr(0,1);
  for(int i=0; i<_n; i++) 
    if(distr(randomNumberGenerator)<=p) v[i]=distr(randomNumberGenerator);
  return v;
}



Vectorh::Vectorh(const Cvector& x): SparseVector(x.n){
  for(int i=0; i<n; i++) if(x(i)!=0) (*this)[i]=x.array[i];}


Vectorh::Vectorh(const Vectorv& x): SparseVector(x.n){
  for(auto& p:x) (*this)(p.first)=p.second;}


Vectorh::Vectorh(const Vectorl& x): SparseVector(x.n){
  for(auto& p:x) (*this)(p.first)=p.second;}


Vectorh::Vectorh(const Vectorl& x, const class Remap& remap, const bool inverse): SparseVector(x.n){
  if(!inverse) for(auto& p:x) (*this)(remap.forward[p.first])=p.second;
  else for(auto& p:x) (*this)(remap.backward[p.first])=p.second;
}


// ---- I/O ------------------------------------------------------------------------------------------------------


string Vectorh::classname(){return "Vectorh";}


Vectorh::Vectorh(Bifstream& ifs): SparseVector(0){
  ifs.check("Vectorh",0);
  ifs.read(n);
  ifs.read_unordered_map(*this);
}


void Vectorh::serialize(Bofstream& ofs) const{
  ofs.tag("Vectorh",0);
  ofs.write(n);
  ofs.write_unordered_map(*this);
}


void Vectorh::serialize(Rstream& rstream) const{
  rstream<<"Vectorh{"<<Rstream::endl;
  rstream.var("n",n);
  for(auto p:*this){rstream<<"("<<p.first<<","<<p.second<<")"<<Rstream::endl;}
  rstream<<"}"<<Rstream::endl;
}


string Vectorh::str(const Sparse dummy) const{ // 
  ostringstream res;
  for(auto& it:*this) res<<it.first<<" : "<<it.second<<endl;
  return res.str();
}








