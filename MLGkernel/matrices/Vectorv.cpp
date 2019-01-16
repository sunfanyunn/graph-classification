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
Vectorv::Vectorv(const int _n, const class Random& dummy): SparseVector(_n), sorted(0){
  uniform_real_distribution<FIELD> distr(0,1);
  for(int i=0; i<n; i++) 
    if(distr(randomNumberGenerator)<=dummy.p) push_back(SVpair(i,distr(randomNumberGenerator)));
}
*/


Vectorv Vectorv::Random(const int _n, const FIELD p){
  Vectorv v(_n); 
  uniform_real_distribution<FIELD> distr(0,1);
  for(int i=0; i<_n; i++) 
    if(distr(randomNumberGenerator)<=p) v.push_back(SVpair(i,distr(randomNumberGenerator)));
  return v;
}


// ---- Conversions -----------------------------------------------------------------------------------------------


Vectorv::Vectorv(const Cvector& x): SparseVector(x.n), sorted(true){
  for(int i=0; i<n; i++) if(x(i)!=0) push_back(SVpair(i,x.array[i]));
}

Vectorv::Vectorv(const Vectorl& x): SparseVector(x.n), sorted(true){
  for(auto& p:x) push_back(p);}

Vectorv::Vectorv(const Vectorh& x): SparseVector(x.n), sorted(false){
  for(auto& p:x) push_back(SVpair(p.first,p.second));}

Vectorv::Vectorv(const Vectorv& x, const class Remap& remap, const bool inverse): SparseVector(x.n), sorted(false){
  if (!inverse) for(auto& p:x) push_back(SVpair(remap.forward[p.first],p.second));
  else for(auto& p:x) push_back(SVpair(remap.backward[p.first],p.second));
}


// ---- Operations -----------------------------------------------------------------------------------------------


// ---- I/O ------------------------------------------------------------------------------------------------------



string Vectorv::classname(){return "Vectorv";}


string Vectorv::str(const Sparse dummy) const{
  ostringstream res;
  for(auto& it:*this) res<<it.first<<" : "<<it.second<<endl;
  return res.str();
}


Vectorv::Vectorv(Bifstream& ifs): SparseVector(0), sorted(0){
  ifs.check("Vectorv",0);
  ifs.read(n);
  ifs.read(sorted);
  ifs.read_vector(*this);
}


void Vectorv::serialize(Bofstream& ofs) const{
  ofs.tag("Vectorv",0);
  ofs.write(n);
  ofs.write(sorted);
  ofs.write_vector(*this);
}


void Vectorv::serialize(Rstream& rstream) const{
  rstream<<"Vectorv{"<<Rstream::endl;
  rstream.var("n",n);
  for(auto p:*this){rstream<<"("<<p.first<<","<<p.second<<")"<<Rstream::endl;}
  rstream<<"}"<<Rstream::endl;
}







