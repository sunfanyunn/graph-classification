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



#include"MLGgraph.hpp"

const int infinity = 2147483647;

void MLGgraph::init(){

  assert(labels.size()==n);
  neighbors=vector< vector<int> >(n);
  subgraphs=vector< unordered_set<int> >(n);
  subinstances.resize(n);

  for(int i=0; i<n; i++) {
    for(int j=0; j<n; j++) {
      if(adj(i,j)>0 || i==j) neighbors[i].push_back(j); 
    }
  }

  for(int i=0; i<n; i++) {
    subgraphs[i].insert(i);
  }
  dist=FloydWarshall(adj);
}

// Given an adjacency matrix, returns a pairwise distance matrix
Cmatrix MLGgraph::FloydWarshall(const Cmatrix& A) const{
  Cmatrix d(n,n);

  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      if (A(i,j)>0) d(i,j) = 1.0 / A(i,j);
      else d(i,j) = infinity;
    }
  }

  for(int i=0; i<n; i++) d(i,i)=0;

  for(int k=0; k<n; k++) {
    for(int i=0; i<n; i++) {
      for(int j=0; j<i; j++) {
	      if(d(i,k) + d(k,j) < d(i,j)) {
          d(i,j) = d(i,k) + d(k,j);
          d(j,i) = d(i,j);
        }
      }
    }
  }

  return d;
}



void MLGgraph::grow_subgraphs(const int radius){
  for(int r=0; r<radius; r++) {
    for(int i=0; i<n; i++) {
      unordered_set<int> newsubgraph=subgraphs[i]; 
      for(auto p: subgraphs[i]) {
	      for(auto q: neighbors[p]) {
	        newsubgraph.insert(q);
        }
      }
      subgraphs[i]=move(newsubgraph);
    }
  }
}


void MLGgraph::double_subgraphs(){
  // for each vertice, look at vertices in its subgraph. take the subgraph of v and add them to this subgraph
  for(int i=0; i<n; i++) {
    unordered_set<int> newsubgraph; 
    for(auto p: subgraphs[i]) {
      for(auto q: subgraphs[p]) {
	      newsubgraph.insert(q);
      }
    }
    subgraphs[i]=move(newsubgraph);
  }
}

Cmatrix MLGgraph::subLaplacian(const vector<int>& v, double eta=0) const{
  int m=v.size();
  Cmatrix L(m,m);
  //vector<int> v(m); {int i=0; for(auto p:subgraph) v[i++]=p;}

  for(int a=0; a<m; a++) {
    double t=0;
    for(int b=0; b<m; b++) {
      if(a!=b) {
        L(a, b) = -adj(v[a], v[b]);
        t += adj(v[a],v[b]);
      }
    }
    L(a,a) = eta + t; // add regularizer eta to the diagonal
  }
  return L;
}

void MLGgraph::push_to_linearizer(Linearizer<FLGinstance>& linearizer, double eta){
  assert(labels.size()==n);
  int nfeatures=labels[0].size();

  // for each vertex of the graph, do something...
  for(int i=0; i<n; i++) {
    int nsub=subgraphs[i].size(); // size of subgraph around this vertex
    vector<int> v(nsub); 
    {
      // create a vector of size of node i's subgraph. for each node in the subgraph
      // basically create a vector of neighbors from a set. order arbitrary?
      int si=0;
      for(auto p:subgraphs[i]) v[si++]=p;
    }
    //{cout<<"subgraph [ ";for(auto p: v) cout<<p<<" "; cout<<"]"<<endl;}
    Cmatrix L=subLaplacian(v, eta);
    Cmatrix U(nsub,nfeatures);
    for(int si=0; si<nsub; si++){
      double d=(dist(i,v[si])+1);
      for(int j=0; j<nfeatures; j++) {
	      U(si,j)=labels[v[si]](j)/d;
      }
    }
    //subinstances[i]=linearizer.add_instance(FLGinstance(move(L),move(U)));
    //cout<<"  Pushing subgraph of size "<<nsub<<endl;  
    subinstances[i]=linearizer.add_instance(move(L),move(U));
  }
}



void MLGgraph::pull_features(){
  for(int i=0; i<n; i++) {
    labels[i]=subinstances[i]->linearization.copy();
    //cout<<labels[i]<<endl;
  }
}


void MLGgraph::compute_flg(){
  vector<int> vertices(n);
  for(int i=0; i<n; i++) vertices[i]=i;
  flg.L=subLaplacian(vertices);

  assert(labels.size()==n);
  int nfeatures=labels[0].size();
  flg.U=Cmatrix(n,nfeatures);
  for(int i=0; i<n; i++) {
    for(int j=0; j<nfeatures; j++) {
	    flg.U(i,j)=labels[i](j);
    }
  }
}


void MLGgraph::computeDegreeFeatures(const int maxdeg){
  for(int j=0; j<n; j++) {
    int d=adj.vcolumn(j).sum(); 
    labels[j]=Cvector::Zero(maxdeg+1);
    if(d<=maxdeg) labels[j](d)=1;
    else labels[j](maxdeg) = 1; // d > maxDeg, and we don't want to have super long feature vectors
  }   
}


string MLGgraph::str() const{
  ostringstream oss;
  oss<<adj.str()<<endl;
  for(int i=0; i<n; i++)
    oss<<labels[i]<<endl;
  return oss.str();
}



    //vector<Cvector> v(subgraphs[i].size());
    //int j=0; for(auto p:subgraphs[i]) v[j++]=(labels[p].copy());
    //subinstances[i]=linearizer.add_instance(FLGinstance(move(L),move(v)));
/*
  int maxdeg=0; 
  for(int j=0; j<n; j++){
    int d=adj.vcolumn(j).sum(); 
    if(d>maxdeg) maxdeg=d;
  }
*/
