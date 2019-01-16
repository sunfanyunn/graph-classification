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

#ifndef _Rstream
#define _Rstream

#include "pMMFbase.hpp"


class Rstream{
public:

  Rstream(ostream& _out, const int _depth=16):out(_out),indent(0),depth(_depth),bol(true){}

  ~Rstream(){out<<std::endl;}

  //const Rstream& indent(const int n=0) const {for(int i=0; i<n; i++) out<<"  "; return *this;}

  template<class T>
  Rstream& operator<<(const T& x){
    //if(typeid(x)==typeid(Rstream::end)) {out<<std::endl; bol=true; return *this;}
    //if(typeid(x)==typeid(endme)) {out<<std::endl; bol=true; return *this;}
    if(bol) {for(int i=0; i<indent; i++) out<<"  "; bol=false;}
    out<<x; return *this;}


  typedef  Rstream& (*RstreamManipulator)(Rstream&);
  Rstream& operator<<(const RstreamManipulator& manip){return manip(*this);}


  static Rstream& endl(Rstream& x){
    x.out<<std::endl; 
    x.bol=true; 
    return x;
  }


  template<class T>
  Rstream& write(const T& x){
    if(depth<0){out<<std::endl; bol=true; return *this;}
    indent++; depth--;
    x.serialize(*this);
    indent--; depth++;
    return *this;
  }


  template<class T>
  const Rstream& var(const char* name, const T& x){
    if(bol) {for(int i=0; i<indent; i++) out<<"  "; bol=false;}
    out<<"  "<<name<<"="<<x<<std::endl; 
    bol=true;
    return *this;
  }



  //typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
  //typedef CoutType& (*StandardEndLine)(CoutType&);
  //Rstream& operator<<(StandardEndLine manip){return *this;}


public:

  int indent;
  int depth;
  mutable bool bol;

  ostream& out;
  
  
};

#endif
