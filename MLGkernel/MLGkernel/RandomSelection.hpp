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



#ifndef _RandomSelection
#define _RandomSelection

#include <unordered_set>


#include "Activemap.hpp"
#include "pMMFbase.hpp"

extern default_random_engine randomNumberGenerator;


class RandomSelection: public vector<int>{
public:

  RandomSelection(const int k, const int n): vector<int>(k){
    assert(k<=n);

    if(k<0.3*n){
      uniform_int_distribution<int> distri(0,n-1);
      for(int i=0; i<k; i++){
	int x; while(selected.find(x=distri(randomNumberGenerator))!=selected.end()){}
	(*this)[i]=x;
	selected.insert(x);
      }
      return;
    }
    
    Activemap amap(n);
    for(int i=0; i<k; i++){
      uniform_int_distribution<int> distri(0,n-i);
      int j=amap(distri(randomNumberGenerator));
      amap.remove(j);
      (*this)[i]=j;
    }
      
  }

public:

  unordered_set<int> selected;
  //Activemap activemap;

};

#endif


      //do{x=distri(randomNumberGenerator);
      //}while(selected.find(x)!=selected.end());
