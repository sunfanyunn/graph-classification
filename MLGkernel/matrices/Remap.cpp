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


#include "Remap.hpp"
#include <random>

extern default_random_engine randomNumberGenerator;


Remap::Remap(const Remap& x): n(x.n){
  forward=new int[n]; for(int i=0; i<n; i++) forward[i]=x.forward[i];
  backward=new int[n]; for(int i=0; i<n; i++) backward[i]=x.backward[i];
}


Remap::Remap(Remap&& x): n(x.n){
  forward=x.forward; x.forward=nullptr;
  backward=x.backward; x.backward=nullptr;
  x.n=0;
}


Remap& Remap::operator=(const Remap& x){
  delete forward; delete backward; n=x.n; 
  forward=new int[n]; for(int i=0; i<n; i++) forward[i]=x.forward[i];
  backward=new int[n]; for(int i=0; i<n; i++) backward[i]=x.backward[i];
  return *this;
}


Remap Remap::operator=(Remap&& x){
  delete forward; delete backward;
  forward=x.forward; x.forward=nullptr;
  backward=x.backward; x.backward=nullptr;
  n=x.n; x.n=0;
  return *this;
}


Remap Remap::Random(const int n){
  Remap R(n);
  for(int i=0; i<n-1; i++){
    uniform_int_distribution<int> distr(i+1,n-1);
    int j=distr(randomNumberGenerator);
    R.swap(i,j);}
  return R;
}


/* DEPRECATED 
Remap::Remap(const int _n, const Random random): Remap(_n){
  for(int i=0; i<n-1; i++){
    uniform_int_distribution<int> distr(i+1,n-1);
    int j=distr(randomNumberGenerator);
    swap(i,j);
  }
}
*/


string Remap::str() const{
  ostringstream stream;
  for(int i=0; i<n; i++) stream<<i<<" -> "<<forward[i]<<endl;
  return stream.str();
}
