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



#include"FLGkernel.hpp"


double FLGkernel::operator()(const FLGinstance& x1, const FLGinstance& x2) const
{

  if(x1.Sinv.nrows==0) const_cast<FLGinstance&>(x1).precompute(gamma);
  if(x2.Sinv.nrows==0) const_cast<FLGinstance&>(x2).precompute(gamma);

  Cvector lambda=(x1.Sinv+x2.Sinv).eigenvalues();
  //double detS=1; for(int i=0; i<lambda.n; i++) detS*=lambda(i); detS=1.0/detS; 
  double log_detS=0; for(int i=0; i<lambda.n; i++) log_detS-=log(lambda(i)); 

  //double r=sqrt(detS/sqrt(x1.detS*x2.detS));
  double logr=(log_detS-0.5*(x1.log_detS+x2.log_detS))/2;
  double r=0;
  if(logr<-30){cout<<"Underflow!"<<endl;} 
  else r=exp(logr);
  /**
  cout << "x1.sinv" << endl;
  cout << x1.Sinv << endl;

  cout << "x2.sinv" << endl;
  cout << x2.Sinv << endl;

  cout << "x1.logdets" << endl;
  cout << x1.log_detS << endl;

  cout << "x2.logdets" << endl;
  cout << x2.log_detS << endl;
  cout << "log dets" << endl;
  cout <<  log_detS << endl;
  **/
  //cout<<"k(.,.)="<<r<<endl;
  return r;
}

