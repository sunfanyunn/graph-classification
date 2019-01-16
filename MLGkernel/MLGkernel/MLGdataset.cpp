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



#include"MLGdataset.hpp"

#include "FLGkernel.hpp"
#include "MatrixOF_ASCII.hpp"
#include "params.hpp"
#include <string>

void MLGdataset::condense(const int nlevels, const int leaf_radius){
  assert(nlevels>0);
  radius = leaf_radius;
  for(int s=0; s<nlevels; s++){
    
    FLGkernel kernel(gamma);
    Linearizer<FLGinstance> linearizer(gamma);

    cout<<"\n---- LEVEL "<<s<<" --------------------------------------------------------------------\n\n";
    cout<<"Growing..."<<endl;
    if(s==0 || (grow == Params::GROW)) {
      for(auto p:graphs) p->grow_subgraphs(leaf_radius);
    } else {
      for(auto p:graphs) p->double_subgraphs();
    }

    cout<<"Extracting..."<<endl;
    for(auto p:graphs) p->push_to_linearizer(linearizer, eta);

    cout<<"Linearizing..."<<endl;
    linearizer.linearize(kernel,50,10);

    cout<<"Pulling features..."<<endl;
    for(auto p:graphs) p->pull_features();
  }

  cout<<"----------------------------------------------------------------------------------\n";
  cout<<"Computing top level instances" <<endl;
  for(auto p:graphs) p->compute_flg();
}


void MLGdataset::computeGram(int levels, int radius){
  condense(levels, radius); 
  int N=graphs.size();
  gram=Cmatrix(N,N);
  FLGkernel kernel(gamma);

  for(int i=0; i<graphs.size(); i++){
    ((*graphs[i]).flg).precompute(gamma);
  }
  cout<<"Computing Gram matrix..."<<endl;
  { int nthreads=threadManager.maxthreads;
    ThreadBank K_threads(nthreads);
    for(int t=0; t<nthreads; t++)
      K_threads.add([this, &kernel, N, nthreads](int t){
		    //{CoutLock lock; cout<<"  Starting thread  "<<t<<endl;}
		    for(int i=t; i<N; i+=nthreads) {
			    for(int j=0; j<=i; j++) {
			      gram(j,i)=(gram(i,j)=kernel(*graphs[i],*graphs[j]));
          }
        }
		    //{CoutLock lock; cout<<"  Finishing thread "<<t<<endl;}
		  },t);
  }

  //cout<<"\nK=\n"<<gram<<endl;
  //cout<<"\nK has "<<gram.nrows << " and " << gram.ncols<<endl;
}


void MLGdataset::loadGraphs(std::string filename){
  ifstream ifs(filename);
  if(ifs.fail()){cout<<"Failed to open "<<filename<<"."<<endl; exit(0);}

  int numGraphs = 0;
  int i=0;
  int n;
  ifs >> numGraphs;
  while(ifs.good()){
    ifs>>n;
    if(!ifs.good()) break;
    //cout<<"Reading graph "<<++i<<" (n="<<n<<")"<<endl;
    Cmatrix A(n,n);
    for(int j=0; j<n; j++){
      for(int k=0; k<n; k++)
	      ifs>>A(j,k);
    }
    graphs.push_back(new MLGgraph(move(A)));
  }
  assert(numGraphs == graphs.size());
}

void MLGdataset::loadDiscreteFeatures(std::string filename, int numFeatures){
  ifstream ifs(filename);
  if(ifs.fail()){
    cout << "Failed to open " << filename << "." << endl;
    exit(0);
  }
  int numVertices = 0;
  int graphIndex = 0;
  int label;
  int numGraphs;
  ifs >> numGraphs;
  while(ifs.good()){
    ifs>>numVertices;
    if(!ifs.good()) break;
    for(int i=0; i<numVertices; i++){
      ifs >> label;
      graphs[graphIndex]->labels[i] = Cvector::Zero(numFeatures+1);
      graphs[graphIndex]->labels[i](label) = 1;
    }
    graphIndex++;
  }
  assert(graphs.size() == numGraphs);
  assert(graphIndex == graphs.size());
}

void MLGdataset::loadFeatures(std::string filename){
// Add features to each graph from a feature ascii file
  ifstream ifs(filename);
  if(ifs.fail()){
    cout<<"Failed to open "<<filename<<"."<<endl;
    exit(0);
  }

  int numGraphs;
  int numVertices = 0;
  int numFeatures = 0;
  int graphIndex = 0;
  ifs >> numGraphs;
  while(ifs.good()){
    ifs>>numVertices;
    ifs>>numFeatures;
    //cout << "Found n: " << numVertices << " and num features: " << numFeatures << endl;
    if(!ifs.good()) break;
    for(int i=0; i<numVertices; i++){
      graphs[graphIndex]->labels[i] = Cvector::Zero(numFeatures);
      for(int j=0; j<numFeatures; j++){
        ifs >> graphs[graphIndex]->labels[i](j);
      }
      //cout << "LABELS FOR GRAPH " << g << " vertex " << i << endl;
      //cout << graphs[g]->labels[i]<<endl;
    }
    graphIndex++;
  }
  assert(graphIndex == graphs.size());
}

void MLGdataset::saveGram(std::string filename){
  cout << "Saving gram to " << filename << endl;
  MatrixOF_ASCII::Dense file(filename,gram.nrows,gram.ncols);
  gram.saveTo(file);
}

/* fillGram: Only used in the python interface. Call it with a numpy matrix
 * of zeros. The given np matrix is then filled with the values from
 * the gram matrix. Assumes that computeGram has been called beforehand
 * and that the size of the numpy matrix is equal to n-by-n, where n = number
 * of graphs in the dataset.
 */
void MLGdataset::fillGram(double *npmatrix, int rows, int cols) {
  // rows == cols == graphs.size()
  int index;
  for(int i=0; i<rows; i++) {
    for(int j=0; j<cols; j++) {
      index = (i * cols) + j;
      npmatrix[index] = gram(i, j);
    }
  }
}
