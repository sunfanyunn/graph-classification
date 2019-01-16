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

#ifndef _Graph
#define _Graph

#include "Cmatrix.hpp"
#include "TopkList.hpp"

extern default_random_engine randomNumberGenerator;


template<class MATRIX>
class Graph{
public:

  Graph(): n(0){}
  Graph(const Graph& G): n(G.n), adj(G.adj){}
  Graph(Graph&& G): n(G.n), adj(move(G.adj)){G.n=0;}
  Graph& operator=(const Graph& G){n=G.n; adj=G.adj; return *this;}
  Graph& operator=(Graph&& G){n=G.n; adj=move(G.adj); G.n=0; return *this;}

  Graph(const int _n):n(_n), adj(n,n){}
  Graph(const int _n, int dummy):n(_n), adj(MATRIX::Zero(n,n)){}

  static Graph<MATRIX> Zero(const int _n);
  static Graph<MATRIX> Path(const int _n);
  static Graph<MATRIX> Cycle(const int _n);
  static Graph<MATRIX> Tree(const int k, const int d);
  static Graph<MATRIX> RegTree(const int k, const int d);
  static Graph<MATRIX> Complete(const int _n);
  static Graph<MATRIX> ErdosRenyi(const int _n, const double p);
  static Graph<MATRIX> Kronecker(const int _d, const double p1, const double p2);
  static Graph<MATRIX> kNN(const int _n, const int k);
  static Graph<MATRIX> Grid(const int m1, const int m2);

public:

  void addEdge(const int i, const int j){adj(i,j)=1; adj(j,i)=1;}

  MATRIX laplacian() const; 

public:

  int n;
  MATRIX adj;


};


template<class MATRIX>
Graph<MATRIX> Graph<MATRIX>::Zero(const int n){
  Graph G(n,0);
  return G;
}


template<class MATRIX>
Graph<MATRIX> Graph<MATRIX>::Path(const int n){
  Graph G(n,0);
  for(int i=0; i<n-1; i++){G.adj(i,i+1)=1; G.adj(i+1,i)=1;}
  return G;
}


template<class MATRIX>
Graph<MATRIX> Graph<MATRIX>::Cycle(const int n){
  Graph G(n,0);
  for(int i=0; i<n; i++){G.adj(i,(i+1)%n)=1; G.adj((i+1)%n,i)=1;}
  return G;
}


template<class MATRIX>
Graph<MATRIX> Graph<MATRIX>::Tree(const int k, const int d){
  if(d==0) return Graph(1,0);
  Graph sub=Tree(k,d-1);
  Graph G(k*sub.n+1,0);
  int offset=1;
  for(int u=0; u<k; u++){
    G.addEdge(0,offset);
    for(int i=0; i<sub.n; i++)
      for(int j=0; j<sub.n; j++)
	G.adj(offset+i,offset+j)=sub.adj(i,j);
    offset+=sub.n;
  }
  return G;
}


template<class MATRIX>
Graph<MATRIX> Graph<MATRIX>::RegTree(const int k, const int d){
  assert(d>=1);
  Graph sub=Tree(k,d-1);
  Graph G((k+1)*sub.n+1,0);
  int offset=1;
  for(int u=0; u<k+1; u++){
    G.addEdge(0,offset);
    for(int i=0; i<sub.n; i++)
      for(int j=0; j<sub.n; j++)
	G.adj(offset+i,offset+j)=sub.adj(i,j);
    offset+=sub.n;
  }
  return G;
  
}

/*
  assert(d>=1);
  Graph G(1+(k+1)*((pow(k,d)-1)/(k-1)),0);
  cout<<G.n<<endl;
  for(int i=0; i<k+1; i++) G.addEdge(0,i+1); 
  int offset=1;
  for(int l=1; l<d; l++){
    int nnodes=(k+1)*pow(k,l-1);
    for(int i=0; i<nnodes; i++)
      for(int u=0; u<k; u++)
	G.addEdge(offset+i,offset+nnodes+k*i+u);
    offset+=nnodes;
  } 
  return G;
}
*/

template<class MATRIX>
Graph<MATRIX> Graph<MATRIX>::Complete(const int n){
  Graph G(n,0);
  for(int i=0; i<n; i++)
    for(int j=0; j<n; j++) 
      if(i!=j) G.adj(i,j)=1;
  return G;
}


template<class MATRIX>
Graph<MATRIX> Graph<MATRIX>::ErdosRenyi(const int n, const double p){
  Graph G(n,0);
  for(int i=0; i<n; i++)
    for(int j=0; j<i; j++){
      uniform_real_distribution<FIELD> distr(0.0,1.0);
      if(distr(randomNumberGenerator)<p) {G.adj(i,j)=1; G.adj(j,i)=1;}
    }
  return G;
}


template<class MATRIX>
Graph<MATRIX> Graph<MATRIX>::Kronecker(const int d, const double p1, const double p2){
  Cmatrix seed(2,2); seed(0,0)=p1; seed(1,1)=p1; seed(0,1)=p2; seed(1,0)=p2;
  Cmatrix P=Cmatrix::Kronecker(seed,d);
  int n=P.nrows;
  Graph G(n,0);
  for(int i=0; i<n; i++)
    for(int j=0; j<i; j++){
      uniform_real_distribution<FIELD> distr(0.0,1.0);
      if(distr(randomNumberGenerator)<P(i,j)) {G.adj(i,j)=1; G.adj(j,i)=1;}
    }
  return G;
}


template<class MATRIX>
Graph<MATRIX> Graph<MATRIX>::Grid(const int m1, const int m2){

  Graph G(m1*m2,0);

  for(int i=0; i<m1-1; i++)
    for(int j=0; j<m2; j++){
      G.adj(i*m2+j,(i+1)*m2+j)=1;
      G.adj((i+1)*m2+j,i*m2+j)=1;
    }
	
  for(int i=0; i<m1; i++)
    for(int j=0; j<m2-1; j++){
      G.adj(i*m2+j,i*m2+j+1)=1;
      G.adj(i*m2+j+1,i*m2+j)=1;
    }
	
  return G;
}


template<class MATRIX>
Graph<MATRIX> Graph<MATRIX>::kNN(const int n, const int k){
  Graph G(n,0);
  Cmatrix X=Cmatrix::Random(n,2);
  for(int i=0; i<n; i++){
    TopkList topk(k);
    for(int j=0; j<n; j++)
      if(j!=i) topk.consider(j,-(pow(X(j,0)-X(i,0),2)+pow(X(j,1)-X(i,1),2)));
    for(auto& p:topk) G.adj(i,p.first)=1;
  }
  for(int i=0; i<n; i++)
    for(int j=0; j<i; j++)
      if(G.adj(i,j)==1||G.adj(j,i)==1){G.adj(i,j)=1; G.adj(j,i)=1;}
  return G;
}


template<class MATRIX>
MATRIX Graph<MATRIX>::laplacian() const{
  MATRIX L(n,n);
  for(int i=0; i<n; i++){
    double t=0; for(int j=0; j<n; j++) if(i!=j) t+=(L(i,j)=-adj(i,j));
    L(i,i)=-t;
  }
  return L;
}


#endif
