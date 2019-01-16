/*
 -----------------------------------------------------------------------------
 
 MLGkernel is an open source implementation of the Multiscale Laplacian Graph
 Kernel for computing the gram matrix of a collection of graphs.
 
 Copyright (C) 2016 Imre Risi Kondor, Horace Pan
 
 
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, see <http://www.gnu.org/licenses/>.
 
 ----------------------------------------------------------------------------- */



#include"FLGinstance.hpp"


void FLGinstance::precompute(const double gamma){
  if(Sinv.nrows>0) return;

  //cout<<"L="<<endl<<L<<endl;
  //cout<<"U="<<endl<<U<<endl;

  Cmatrix Linv=invert(L);
  //cout<<"Linv="<<endl<<Linv<<endl;

  //int n=Linv.nrows;
  //assert(labels.size()==n);
  //int m=labels[0].n;

  //Cmatrix U(n,m);
  //for(int i=0; i<n; i++)
  //  for(int j=0; j<m; j++)
  //    U(i,j)=labels[i](j);
  Cmatrix S=U.dot(Linv*U);
  for(int i=0; i<S.nrows; i++) S(i,i)+=gamma;
  //cout<<"S="<<endl<<S<<endl;

  //Sinv=invert(S,0,&detS)*0.5;
  Sinv=invert(S,0,&log_detS)*0.5;
  //cout<<"Sinv="<<endl<<Sinv<<endl;
}



Cmatrix FLGinstance::invert(const Cmatrix& M, const int _maxrank, double* detp) const{

  pair<Cmatrix*,Cvector*> eigenp=M.symmetricEigensolver();
  Cmatrix& eigs=*eigenp.first;
  Cvector& lambda=*eigenp.second;
  int n=M.nrows;
  //cout<<"eigs="<<eigs<<endl;

  int maxrank=min(_maxrank,n);
  if(maxrank==0) maxrank=n;
  for(int i=0; i<lambda.n-1; i++)
    if(lambda(i+1)<lambda(i)) cout<<"WARNING: eigenvalues not sorted in FLGinstance::invert"<<endl;
  for(int i=0; i<maxrank; i++)
    if(lambda(n-1-i)<10e-5*lambda(n-1)) maxrank=i;
  if(detp != nullptr) maxrank = n;
  Cmatrix R=Cmatrix::Zero(n,n);
  //cout<<"maxrank="<<maxrank<<endl;
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++)
      for(int p=0; p<maxrank; p++)
	R(i,j)+=eigs(i,n-1-p)*eigs(j,n-1-p)/lambda(n-1-p);  // CORRECTED

  //if(detp!=nullptr){
  //  *detp=1.0; for(int i=0; i<lambda.n; i++) (*detp)*=lambda(i);
  //}
  if(detp!=nullptr){
    *detp=0.0; for(int i=0; i<lambda.n; i++) (*detp)+=log(lambda(i));
  }

  delete eigenp.first;
  delete eigenp.second;

  return R;
}


