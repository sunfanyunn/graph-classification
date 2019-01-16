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

#ifndef _Log
#define _Log

#include <chrono>
#include "pMMFbase.hpp"

class LogStream{
public:

  //virtual void write(const string& s)=0;
  virtual LogStream& operator<<(const char* s)=0;
  virtual LogStream& operator<<(const string& s)=0;
  virtual LogStream& operator<<(const int& x)=0;
  virtual LogStream& operator<<(const double& x)=0;

};


class Log{
public:

  Log(){startClock();}

public:

  Log& operator<<(const string& s);
  Log& operator<<(const char* s);

  /*
  Log& skipline(const int n=1, const int v=0){
    if(verbosity<v) return *this;
    if(skippedlines>=n) return *this;
    for(int i=0; i<n-skippedlines; i++) cout<<endl; skippedlines=n;
    return *this;
  }
  */

  Log& skip(const int v=0, const int n=1){
    if(verbosity<v) return *this;
    if(skippedlines>=n) return *this;
    if(stream==nullptr) for(int i=0; i<n-skippedlines; i++) cout<<endl; 
    else for(int i=0; i<n-skippedlines; i++) (*stream)<<"\n";
    skippedlines=n;
    return *this;
  }

  Log& log(const int v, const string& s);
  Log& log(const int v, const char* s, const int i);
  Log& log(const int v, const char* s1, const int i, const char* s2);
  Log& log(const int v, const char* s1, const int i1, const char* s2, const int i2, const char* s3);
  Log& log(const int v, const char* s1, const int i1, const char* s2, const int i2, const char* s3, const int i3, const char* s4);
  Log& log(const int v, const char* s1, const int i1, const char* s2, const int i2, const char* s3, const int i3, const char* s4, const int i4, const char* s5);
  Log& log(const int v, const char* s, const double f);
  Log& log(const int v, const int i, const double f);

  void startClock(const int i=0);
  double clock(const int i=0);

public:
  
  //vector<chrono::time_point<chrono::system_clock> > time;
  
  chrono::time_point<chrono::system_clock> t;

  int verbosity=0;
  int skippedlines=0;

  LogStream* stream=nullptr;
};



#endif
