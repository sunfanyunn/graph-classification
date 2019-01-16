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



#ifndef _MLGgraph
#define _MLGgraph

#include <unordered_set>

#include "Cmatrix.hpp"
#include "Graph.hpp"
#include "FLGinstance.hpp"
#include "Linearizer.hpp"


class MLGgraph{
public:

  MLGgraph(const MLGgraph& x): n(x.n), adj(x.adj.copy()), labels(x.labels.size()){
    for(int i=0; i<x.labels.size(); i++) labels[i]=x.labels[i].copy(); init();}
  

public:

  MLGgraph(Cmatrix&& _adj): adj(move(_adj)){
    n=adj.nrows; labels=vector<Cvector>(n); for(auto& p:labels) p=Cvector::Filled(1,0); init();
  }
  MLGgraph(Graph<Cmatrix>&& G){
    n=G.n; adj=move(G.adj); labels=vector<Cvector>(n); init();
  }
  MLGgraph& operator=(Graph<Cmatrix>&& G){
    n=G.n; adj=move(G.adj); labels=vector<Cvector>(n); init(); return *this;
  }
  
public:

  void grow_subgraphs(const int radius);
  void double_subgraphs();
  void push_to_linearizer(Linearizer<FLGinstance>& linearizer, double eta);
  void pull_features();
  void compute_flg();

  void computeDegreeFeatures(const int maxdeg);

  string str() const;

private:
  
  void init();
  Cmatrix subLaplacian(const vector<int>& vset, double eta) const;
  Cmatrix FloydWarshall(const Cmatrix& A) const;

public:

  int n;
  Cmatrix adj;
  vector<Cvector> labels;

  vector< vector<int> > neighbors;
  vector< unordered_set<int> > subgraphs;
  vector<Lwrapper<FLGinstance>*> subinstances;

  FLGinstance flg;
  
  Cmatrix dist;
  

};


namespace std{
template<>
class hash< Hwrapper<FLGinstance> >{
public:
  size_t operator()(const Hwrapper<FLGinstance>& x) const{
    return hash<FLGinstance>()(*x.ptr);}
};
};


#endif
