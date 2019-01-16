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
#ifndef _MatrixIF_ASCII
#define _MatrixIF_ASCII

#include "MatrixIF.hpp"


class MatrixIF_ASCII: public MatrixIF{
public:
  class Dense;
  class Sparse;

};


class MatrixIF_ASCII::Dense: public MatrixIF_ASCII{
public:

  Dense(const string filename){
    sparse=0;
    ifs.open(filename);
    if(ifs.fail()){cout<<"Failed to open "<<filename<<"."<<endl; return;}
    char buffer[1024]; 
    int linelength=0; 
    do{
      ifs.get(buffer,1024);
      linelength+=strlen(buffer);
    }while(strlen(buffer)>0);
    // cout<<"Line length="<<linelength<<endl;
    ifs.close(); ifs.open(filename); // why does ifs.seekg(0) not work?
    float b;
    ncols=0; while(ifs.good() && ifs.tellg()<linelength){ifs>>b; ncols++;}
    nrows=0; while(ifs.good()) {for(int i=0; i<ncols;i++) ifs>>b; nrows++;}
    // cout<<nrows<<" "<<ncols<<endl;
    ifs.close();
    ifs.open(filename); i=0; j=0;
  }


public:

  void rewind(){i=0; j=0; }

   MatrixIF& operator>>(IndexValueTriple& dest){
   dest.i=i; dest.j=j;
   if(++j>=ncols) { j=0; i++; }
   if(ifs.good() && i<=nrows) ifs>>dest.value; else {dest.i=-1; return *this;}
   return *this;
  }

  MatrixIF& operator>>(FIELD& dest){
    if(++j>=ncols) {  j=0; i++; }
    if(ifs.good() && i<=nrows) ifs>>dest;
    return *this;
  }
  

public:

  int i;
  int j;
  bool eof;

};



class MatrixIF_ASCII::Sparse: public MatrixIF_ASCII{
public:

  Sparse(const string filename){
    sparse=1;
    ifs.open(filename);
    char buffer[255];
    ifs.get(buffer,255);
    ifs.close(); 

    ifs.open(filename);
    int nextracted=0;
    while(ifs.good() && ifs.tellg()<strlen(buffer)){float b; ifs>>b; if(!ifs.fail()) nextracted++;}
    if(nextracted==2){ifs.close(); ifs.open(filename); ifs>>nrows>>ncols; return;}
    if(nextracted==3){
      ifs.close(); ifs.open(filename);
      nrows=0; ncols=0;
      int a; int b; float f;
      while(ifs.good()){
	ifs>>a>>b>>f;
	if(a>nrows-1) nrows=a+1;
	if(b>ncols-1) ncols=b+1;
      }
      ifs.close(); ifs.open(filename);
      return;
    }
    cout<<"Error: could not parse first line"<<endl;
  }


public:

  MatrixIF& operator>>(IndexValueTriple& dest){
    if(!ifs.good()){dest.i=-1; return *this;}
    ifs>>dest.i>>dest.j>>dest.value; return *this;
  }

public:

};


#endif
