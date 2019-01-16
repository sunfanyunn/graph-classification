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


#include "Cvector.hpp"
#include "Vectorv.hpp"
#include "Vectorl.hpp"
#include "Vectorh.hpp"

#ifdef _withEigen
#include "EigenInterface.hpp"
#endif

extern default_random_engine randomNumberGenerator;


/*
Cvector::Cvector(const int _n, const class Random& dummy): DenseVector(_n) {
  array=new FIELD[n]; 
  uniform_real_distribution<FIELD> distr;
  for(int i=0; i<n; i++) array[i]=distr(randomNumberGenerator);
}
*/

Cvector Cvector::Random(const int n){
  Cvector v(n);
  uniform_real_distribution<FIELD> distr;
  for(int i=0; i<n; i++) v.array[i]=distr(randomNumberGenerator);
  return v;
}

Cvector::Cvector(const initializer_list<FIELD> list): DenseVector(list.size()){
  array=new FIELD[n]; int i=0; for(FIELD v:list) array[i++]=v;
}


Cvector::Cvector(const int _n, const FIELD* _array): DenseVector(_n){
  array=new FIELD[n]; for(int i=0; i<n; i++) array[i]=_array[i];
}



// ---- Conversions ---------------------------------------------------------------------------------------------


/*
Cvector::Cvector(const Vectorv& x): DenseVector(x.n){
  array=new FIELD[n]; for(int i=0; i<n; i++) array[i]=0;
  for(auto& it:x) array[it.first]=it.second;
}


Cvector::Cvector(const Vectorl& x): DenseVector(x.n){
  array=new FIELD[n]; for(int i=0; i<n; i++) array[i]=0;
  for(auto& it:x) array[it.first]=it.second;
}


Cvector::Cvector(const Vectorh& x): DenseVector(x.n){
  array=new FIELD[n]; for(int i=0; i<n; i++) array[i]=0;
  for(auto& it:x) array[it.first]=it.second;
}
*/

#ifdef _withEigen
Cvector::Cvector(const EigenVectorXdAdaptor& x):Cvector(x.size()){
  for(int i=0; i<n; i++) array[i]=x(i);
}
#endif


#ifdef _withEigen
template<>
Eigen::VectorXd Cvector::convert() const{
  Eigen::VectorXd v(n);
  for(int i=0; i<n; i++) v(i)=array[i];
  return v;
}
#endif


Cvector Cvector::merge(const Cvector& x, const Cvector& y){
  Cvector r(x.n+y.n);
  for(int i=0; i<x.n; i++) r.array[i]=x.array[i];
  for(int i=0; i<y.n; i++) r.array[i+x.n]=y.array[i];
  return r;
}

Cvector::Virtual Cvector::vsubvector(const int beg, const int end){
  return Cvector::Virtual(end-beg,&array[beg]);}


// ---- I/O ------------------------------------------------------------------------------------------------------



ostream& operator<<(ostream& stream, const Cvector& x){
  stream<<"("; for (int i=0; i<x.n-1; i++) stream<<x.array[i]<<","; stream<<x.array[x.n-1]<<")";
  return stream;
}


string Cvector::classname(){return "Cvector";}


void Cvector::serialize(Rstream& rstream) const{
  rstream<<"Cvector{"<<Rstream::endl;
  rstream.var("n",n);
  rstream.var("array"," *OMITTED*");
  rstream<<"}"<<Rstream::endl;
}


void Cvector::serialize(Bofstream& ofs) const{
  ofs.tag("Cvector",0);
  ofs.write(n);
  ofs.write_array(array,n);
}


Cvector::Cvector(Bifstream& ifs): DenseVector(0){
  ifs.check("Cvector",0);
  ifs.read(n);
  ifs.read_array(array);
}
