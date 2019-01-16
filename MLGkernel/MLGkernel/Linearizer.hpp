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



#ifndef _Linearizer
#define _Linearizer

#include "RandomSelection.hpp"
#include "Activemap.hpp"
#include "Cmatrix.hpp"
#include "FLGinstance.hpp"
#include "ThreadBank.hpp"

extern ThreadManager threadManager;


template<class INSTANCE>
class Lwrapper{
public:
  Lwrapper(const INSTANCE* _ptr): ptr(_ptr){}
  ~Lwrapper(){
    delete ptr;
  }
  const INSTANCE* ptr;
  Cvector linearization;
};


template<class INSTANCE>
class Hwrapper{
public:
  Hwrapper(const INSTANCE* _ptr): ptr(_ptr){}
  bool operator==(const Hwrapper& x) const {return (*ptr)==*(x.ptr);}
  const INSTANCE* ptr;
};



template<class INSTANCE>
class Linearizer{
public:
  Linearizer(double g): gamma(g) {};
  ~Linearizer() {
    for(auto p:uniques) delete p;
  }
  
public:

  template<class ARG1, class ARG2>
  Lwrapper<INSTANCE>* add_instance(ARG1&& arg1, ARG2&& arg2){
    INSTANCE* x=new INSTANCE(move(arg1),move(arg2)); // take ownership
    Hwrapper<INSTANCE> hwrap(x);
    auto it=htable.find(hwrap);
    if(it!=htable.end()){
      cache_hits++;
      delete x; 
      return it->second;
    }
    auto* lwrap=new Lwrapper<INSTANCE>(x);
    uniques.push_back(lwrap);
    htable[hwrap]=lwrap;
    return lwrap;
  }

  template<class KERNEL>
  void linearize(const KERNEL& kernel, const int nsample, const int maxdim);
  
  template<class KERNEL>
  void project_instances(const KERNEL& kernel, vector<Lwrapper<INSTANCE>*>& _uniques);

  template<class KERNEL>
  void project(const KERNEL& kernel, Lwrapper<INSTANCE>& x) const;


public:

  vector< Lwrapper<INSTANCE>* > uniques;
  unordered_map< Hwrapper<INSTANCE>, Lwrapper<INSTANCE>* > htable;
  
  int dim=0;
  int nanchors=0;
  vector<int> anchors;
  Cmatrix Q;
  double gamma;
  int cache_hits=0;
};



template<class INSTANCE>
template<class KERNEL>
void Linearizer<INSTANCE>::project(const KERNEL& kernel, Lwrapper<INSTANCE>& x) const{
  Cvector kvec(nanchors);
  for(int i=0; i<nanchors; i++) {
    kvec(i)=kernel(*x.ptr,*uniques[anchors[i]]->ptr);
  }
  // Q is the matrix of eigenvalues of the reduced Kernel
  x.linearization=Q.dot(kvec);
  //cout<<x.linearization<<endl;
}


template<class INSTANCE>
template<class KERNEL>
void Linearizer<INSTANCE>::linearize(const KERNEL& kernel, const int nsample, const int maxdim){
  // Creates the kernel matrix/reduced kernel matrix and then projects things
  //cout<<"Linearizing..."<<endl;

  int num_sample = nsample;
  if(uniques.size() < num_sample) num_sample=uniques.size();
  vector<int> sel=RandomSelection(num_sample,uniques.size());
  // Precompute everything since we need to do this anyway
  for(int i=0; i<uniques.size(); i++){
    if((*uniques[i]->ptr).Sinv.nrows == 0)
      const_cast<FLGinstance&>(*uniques[i]->ptr).precompute(gamma);
  }
  // Fill the kernel matrix K
  cout<<"Computing kernel...\n";
  Cmatrix K(num_sample, num_sample);
  {
    int nthreads=threadManager.maxthreads;
    ThreadBank K_threads(nthreads);
    for(int t=0; t<nthreads; t++)
      K_threads.add([this, &kernel, &K, &sel, num_sample, nthreads](int t){
		    //{CoutLock lock; cout<<"  Starting thread  "<<t<<endl;}
		    for(int i=t; i<num_sample; i+=nthreads) {
			    for(int j=0; j<=i; j++) {
		        //{CoutLock lock; cout<<"thread " << t << " filling entry " << i << ", " << j <<endl;}
			      K(j,i)=(K(i,j)=kernel(*uniques[sel[i]]->ptr,*uniques[sel[j]]->ptr));
          }
        }
		    //{CoutLock lock; cout<<"  Finishing thread "<<t<<endl;}
      },t);
  }
  cout << "Done computing kernel" << endl;
  //cout<<K<<endl;
  //cout<<"  nuniques = "<<uniques.size()<<endl;
  //cout<<"  cachehit = "<<cache_hits<<endl;
  //cout<<"  nsampled = "<<num_sample<<endl;
  double threshold=(K.norm2()/num_sample)*10e-08;
  Activemap amap(num_sample);

  // Remove columns that are too closely correlated with another column
  for(int i=0; i<num_sample; i++) {
    auto coli=K.vcolumn(i);
    for(int j=0; j<i; j++) {
      if(!amap.isactive(j)) continue; 
      auto colj=K.vcolumn(j);
      if((coli-colj).norm2()<threshold) {amap.remove(i); break;}
    }
  }

  //cout<<"  nanchors = "<<amap.nactive<<endl;
  nanchors=amap.nactive; // number sampled - cols that are too correlated
  anchors.resize(nanchors); // anchors is the stuff sampled?? sel is the stuff sampled
  for(int i=0; i<nanchors; i++) anchors[i]=sel[amap(i)]; // vector of the sampled points' id
  Cmatrix Kred(nanchors,nanchors); // reduced kernel matrix after pruning correlated columns
  for(int i=0; i<nanchors; i++)
    for(int j=0; j<nanchors; j++)
      Kred(i,j)=K(amap(i),amap(j));
  //cout<<Kred<<endl;

  pair<Cmatrix*,Cvector*> eigp=Kred.symmetricEigensolver();
  Cmatrix& eigs=*eigp.first;
  Cvector& lambda=*eigp.second;

  for(int i=0; i<lambda.n-1; i++)
    if(lambda(i+1)<lambda(i)) cout<<"WARNING: eigenvalues not sorted in FLGinstance::invert"<<endl;

  // We want the square root of the eigenvalues
  dim=min(maxdim,nanchors);
  for(int i=0; i<dim; i++)
    if(lambda(nanchors-1-i)<10e-3*lambda(nanchors-1)) dim=i;
  for(int i=0; i<dim; i++)
    lambda(nanchors-1-i)=sqrt(lambda(nanchors-1-i)); // lambdas really holds sqrt ofeigvals
  
  //cout<<"  new dim  = "<<dim<<endl;
  //cout<<"  lambda   = "<<lambda.vsubvector(nanchors-dim,nanchors)<<endl;
  // take the top-{dim} eigenvectors of the reduced Kernel matrix.
  Q=Cmatrix(nanchors,dim);
  for(int i=0; i<nanchors; i++)
    for(int j=0; j<dim; j++)
      Q(i,j)=eigs(i,nanchors-1-j)/lambda(nanchors-1-j);
  
  cout<<"Projecting ..."<<endl;
  project_instances(kernel, uniques);
  delete eigp.first;
  delete eigp.second;
}

template<class INSTANCE>
template<class KERNEL>
void Linearizer<INSTANCE>::project_instances(const KERNEL& kernel, vector<Lwrapper<INSTANCE>*>& _uniques){
  int nthreads=threadManager.maxthreads;
  int chunk=_uniques.size()/nthreads;
  int ninst=chunk;
  ThreadBank projector_threads(nthreads);
  for(int t=0; t<nthreads; t++){
    if(t==nthreads-1) ninst=_uniques.size()-t*chunk;
    projector_threads.add([this, &kernel, &_uniques](int t, int chunk, int ninst){
      for(int i=0; i<ninst; i++)
        project(kernel,*_uniques[t*chunk+i]);
    },t,chunk,ninst);
  }
}


#endif
