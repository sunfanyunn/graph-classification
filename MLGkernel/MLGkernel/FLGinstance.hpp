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



#ifndef _FLGinstance
#define _FLGinstance

#include "Cmatrix.hpp"
#include "Cvector.hpp"


class FLGinstance{
public:

  //FLGinstance(Cmatrix&& _L, vector<Cvector>&& _labels): 
  //  L(move(_L)), labels(_labels){};

  FLGinstance(){}

  FLGinstance(Cmatrix&& _L, Cmatrix&& _U): 
    L(move(_L)), U(move(_U)){};


public:

  void precompute(const double gamma);

  bool operator==(const FLGinstance& x) const{
    if(L!=x.L) return false;
    if(U!=x.U) return false;
    //if(labels.size()!=x.labels.size()) return false;
    //if(labels!=x.labels) return false;
    return true;
  }

  string str(){
    ostringstream oss; oss<<L<<U<<endl; return oss.str();}

private:

  Cmatrix invert(const Cmatrix& M, const int _maxrank=0, double* detp=nullptr) const;

public:

  Cmatrix L;
  //vector<Cvector> labels;
  Cmatrix U;

  Cmatrix Sinv; // actually Sinv/2
  //double detS;
  double log_detS;

  // Cvector linearization;

};




namespace std{
  template<>
  class hash<FLGinstance>{
  public:
    size_t operator()(const FLGinstance& x) const{
      size_t h=hash<Cmatrix>()(x.L)^hash<Cmatrix>()(x.U);
      //for(auto& p: G.labels) h=(h<<1)^hash<Cvector>()(p); 
      return h;
    }
  };
};


#endif
