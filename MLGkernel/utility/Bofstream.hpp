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

#ifndef _Bofstream
#define _Bofstream

#include "pMMFbase.hpp"
#include <fstream>
#include <string.h>
#include <list>
#include <unordered_map>


class Bofstream: public ofstream{
public:

  Bofstream(const string filename){
    open(filename,ios::binary);
  }

  ~Bofstream(){close();}

public:

  Bofstream& write(const int& x) {ofstream::write(reinterpret_cast<const char*>(&x),sizeof(int)); return *this;}
  Bofstream& write(const double& x) {ofstream::write(reinterpret_cast<const char*>(&x),sizeof(double)); return *this;}
  Bofstream& write(const bool& x) {ofstream::write(reinterpret_cast<const char*>(&x),sizeof(bool)); return *this;}
  Bofstream& write(const BlockIndexPair& x){
    ofstream::write(reinterpret_cast<const char*>(&x.block),sizeof(INDEX));
    ofstream::write(reinterpret_cast<const char*>(&x.index),sizeof(INDEX)); return *this;}
  Bofstream& write(const SVpair& x){
    ofstream::write(reinterpret_cast<const char*>(&x.first),sizeof(INDEX));
    ofstream::write(reinterpret_cast<const char*>(&x.second),sizeof(FIELD)); return *this;}

  template<class TYPE>
  Bofstream& write_vector(const vector<TYPE>& x){
    int n=x.size();
    ofstream::write(reinterpret_cast<const char*>(&n),sizeof(int)); 
    for(auto p:x) write(p);
    return *this;
  }

  template<class TYPE>
  Bofstream& write_list(const list<TYPE>& x){
    int n=x.size();
    ofstream::write(reinterpret_cast<const char*>(&n),sizeof(int)); 
    for(auto p:x) write(p);
    return *this;
  }

  template<class KEY, class VAL>
  Bofstream& write_unordered_map(const unordered_map<KEY,VAL>& x){
    int n=x.size();
    ofstream::write(reinterpret_cast<const char*>(&n),sizeof(int)); 
    for(auto p:x){
      ofstream::write(reinterpret_cast<const char*>(&p.first),sizeof(KEY));
      ofstream::write(reinterpret_cast<const char*>(&p.second),sizeof(VAL));
    }
    return *this;
  }

  template<class TYPE>
  Bofstream& write_array(const TYPE* x, const int n){
    ofstream::write(reinterpret_cast<const char*>(&n),sizeof(int)); 
    ofstream::write(reinterpret_cast<const char*>(x),n*sizeof(TYPE)); 
    //for(int i=0; i<n; i++) write(x[i]);
    return *this;
  }


public:

  template<class TYPE>
  Bofstream& seriate(const TYPE* x){x->serialize(*this); return *this;}

  template<class TYPE>
  Bofstream& serialize(const TYPE& x){x.serialize(*this); return *this;}

  template<class TYPE>
  Bofstream& seriate_vector(const vector<TYPE*>& x){
    int n=x.size();
    ofstream::write(reinterpret_cast<const char*>(&n),sizeof(int)); 
    for(auto p:x) p->serialize(*this);
    return *this;
  }

 template<class TYPE>
  Bofstream& seriate_vector(const vector<TYPE>& x){
    int n=x.size();
    ofstream::write(reinterpret_cast<const char*>(&n),sizeof(int)); 
    for(auto p:x) p.serialize(*this);
    return *this;
  }

  template<class TYPE>
  Bofstream& seriate_array(TYPE** x, const int n){ // why did const have to be removed?
    ofstream::write(reinterpret_cast<const char*>(&n),sizeof(int)); 
    for(int i=0; i<n; i++) x[i]->serialize(*this);
    return *this;
  }


public:

  void tag(const char* name, const int _version){
    int version=_version;
    ofstream::write(name,strlen(name));
    ofstream::write(reinterpret_cast<const char*>(&version),sizeof(int)); 
  }


};


#endif


/*
  template<class TYPE>
  Bofstream& seriate_vector(const vector<TYPE>& x){
    int n=x.size();
    ofstream::write(reinterpret_cast<const char*>(&n),sizeof(int)); 
    for(auto p:x) p.serialize(*this);
    return *this;
  }

  template<class TYPE>
  Bofstream& seriate_array(const TYPE*& x, const int n){
    ofstream::write(reinterpret_cast<const char*>(&n),sizeof(int)); 
    for(int i=0; i<n; i++) x[i].serialize(*this);
    return *this;
  }

 */

/*
  Bofstream& write(const vector<BlockIndexPair>& x){
    int sizet=x.size();
    ofstream::write(reinterpret_cast<const char*>(&sizet),sizeof(int)); 
    for(auto p:x) write(p);
    return *this;
  }
*/
