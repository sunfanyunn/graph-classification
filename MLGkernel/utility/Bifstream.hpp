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

#ifndef _Bifstream
#define _Bifstream

#include "pMMFbase.hpp"
#include <fstream>
#include <string.h>
#include <list>
#include <unordered_map>


extern char strbuffer[255];


class Bifstream: public ifstream{
public:

  Bifstream(const string filename){
    open(filename,ios::binary);
  }

  ~Bifstream(){close();}

public:

  Bifstream& read(int& x) {ifstream::read(reinterpret_cast<char*>(&x),sizeof(int)); return *this;}
  Bifstream& read(double& x) {ifstream::read(reinterpret_cast<char*>(&x),sizeof(double)); return *this;}
  Bifstream& read(bool& x) {ifstream::read(reinterpret_cast<char*>(&x),sizeof(bool)); return *this;}
  Bifstream& read(BlockIndexPair& x){
    ifstream::read(reinterpret_cast<char*>(&x.block),sizeof(INDEX));
    ifstream::read(reinterpret_cast<char*>(&x.index),sizeof(INDEX)); return *this;}
  Bifstream& read(SVpair& x){
    ifstream::read(reinterpret_cast<char*>(&x.first),sizeof(INDEX));
    ifstream::read(reinterpret_cast<char*>(&x.second),sizeof(FIELD)); return *this;}

  template<class TYPE>
  Bifstream& read_vector(vector<TYPE>& x){
    int n; ifstream::read(reinterpret_cast<char*>(&n),sizeof(int)); x.resize(n);
    for(int i=0; i<n; i++) read(x[i]);
    return *this;
  }

  template<class TYPE>
  Bifstream& read_list(list<TYPE>& x){
    int n; ifstream::read(reinterpret_cast<char*>(&n),sizeof(int)); 
    for(int i=0; i<n; i++) {TYPE z; read(z); x.push_back(z);}
    return *this;
  }

  template<class KEY, class VAL>
  Bifstream& read_unordered_map(unordered_map<KEY,VAL>& x){
    int n; ifstream::read(reinterpret_cast<char*>(&n),sizeof(int)); 
    for(int i=0; i<n; i++){
      KEY a; ifstream::read(reinterpret_cast<char*>(&a),sizeof(KEY));
      VAL b; ifstream::read(reinterpret_cast<char*>(&b),sizeof(VAL));
      x[a]=b;}
    return *this;
  }

  template<class TYPE>
  Bifstream& read_array(TYPE*& x){
    int n; ifstream::read(reinterpret_cast<char*>(&n),sizeof(int)); x=new TYPE[n];
    ifstream::read(reinterpret_cast<char*>(x),n*sizeof(TYPE)); 
    return *this;
  }


public:

  template<class TYPE>
  Bifstream& unseriate(TYPE*& x){x=new TYPE(*this); return *this;}

  template<class TYPE>
  Bifstream& unserialize(TYPE& x){x=TYPE(*this); return *this;}

  template<class TYPE>
  Bifstream& unseriate_vector(vector<TYPE*>& x){
    int n; ifstream::read(reinterpret_cast<char*>(&n),sizeof(int)); x.resize(n);
    for(int i=0; i<n; i++) x[i]=new TYPE(*this);
    return *this;
  }

  template<class TYPE>
  Bifstream& unseriate_vector(vector<TYPE>& x){
    int n; ifstream::read(reinterpret_cast<char*>(&n),sizeof(int)); x.resize(n);
    for(int i=0; i<n; i++) x[i]=static_cast<TYPE>(*this);
    return *this;
  }

  template<class TYPE>
  Bifstream& unseriate_array(TYPE**& x){
    int n; ifstream::read(reinterpret_cast<char*>(&n),sizeof(int)); x=new TYPE*[n];
    for(int i=0; i<n; i++) x[i]=new TYPE(*this);
    return *this;
  }


public:

  bool checkname(char const* name) {
    peek().copy(strbuffer,255,0);
    if(!strcmp(strbuffer,name)) return true;
    return false;
  }

  bool check(char const* name, const int version){

    ifstream::get(strbuffer,255,'\0');
    if(strcmp(strbuffer,name)!=0){
      cout<<"Error in Bifstream: "<<name<<" expected, "<<strbuffer<<" found."<<endl; return true;}

    int readver; ifstream::read(reinterpret_cast<char*>(&readver),sizeof(int));
    if(readver!=version){
      cout<<"Error in Bifstream: version "<<version<<" expected, "<<readver<<" found."<<endl; return true;}

    return false;
  }

  string peek(){
    streampos pos=ifstream::tellg();
    ifstream::get(strbuffer,255,'\0');
    ifstream::seekg(pos);
    return strbuffer;
 }

  

public:

  //  static char buffer[255];


};

#endif
