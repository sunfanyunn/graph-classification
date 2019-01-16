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


#ifndef _pMMFbase 
#define _pMMFbase 


// Preprocessor variables  ------------------------------------------------------------------------


// #define _withEigen defined from command line 

// Comment out the following to suppress warnings when objects are copied or assigned
#define _UTILCOPYWARNING
#define _MATRIXCOPYWARNING
#define _BLOCKEDCOPYWARNING 
#define _MCOPYWARNING
#define _MMFCOPYWARNING
#define _MSSCOPYWARNING

// Comment out the following to suppress warnings when objects are moved or move-assigned
// #define _UTILMOVEWARNING
// #define _MATRIXMOVEWARNING
// #define _BLOCKEDMOVEWARNING 
// #define _MMOVEWARNING 
// #define _MMFMOVEWARNING
// #define _MSSMOVEWARNING

// comment out experimental features below
#define _FASTGRAM


// Typedefs ---------------------------------------------------------------------------------------


typedef int INDEX;
typedef double FIELD;


// STL libraries ----------------------------------------------------------------------------------


#include <iostream>
#include <vector>
#include <sstream>
#include <random>
#include <initializer_list>
#include <assert.h>
#include <string>
#include <cstring>
#include <algorithm>
#include <mutex>
#include <atomic>

using namespace std;


// Forward declarations ----------------------------------------------------------------------------


class Zero{};
class Identity{};
class Uninitialized{};
class Dangling{};
class Nullcols{};
class Dense{};
class Sparse{};
class Symmetric{};
class Transpose{};

class Remap;
class Rotation;
class ElementaryRotation;
class GivensRotation;
class KpointRotation; 
template<class TYPE> class Clusters;

class Vector;
class Cvector;
class Vectorv;
class Vectorl;
class Vectorh;

class Matrix;
class Cmatrix; 
class AtomicCmatrix;
template<class COLUMNTYPE> class MatrixX;

template<class BLOCK> class BlockedVector;
template<class BLOCK> class Tower;
template<class BLOCK> class Street;
template<class BLOCK> class BlockedMatrix;

class SparseMatrixFile; 
class DenseMatrixFile;

class MatrixIF;
class MatrixOF;

class Bofstream;
class Bifstream;

template<class BLOCK> class MMFprocess;
template<class BLOCK> class MMFmatrix;
class MMFchannel;
class MMFstage;
class MMF;

template<class BLOCK> class MSSprocess;
template<class BLOCK> class MSSmatrix;
class MSSstage;
class MSS;


// Helper classes ---------------------------------------------------------------------------------


class Random{
public: 
  Random(){} 
  Random(const double _p):p(_p){} 
  Random(const int _n):n(_n){}
public:
  double p=0.5;
  int n=0;
};

template<class T1, class T2, class T3>
struct triple{
public:
  triple(){}
  triple(const T1& v1, const T2& v2, const T3& v3): first(v1), second(v2), third(v3){}
public:
  T1 first; T2 second; T3 third;
};
  

struct IndexValuePair{
  IndexValuePair():i(0),value(0){}
  string str(){ostringstream result; result<<"("<<i<<","<<value<<")"; return result.str();}
  INDEX i; FIELD value; };

struct IndexValueTriple{
  IndexValueTriple():i(0),j(0),value(0){}
  IndexValueTriple(const int _i, const int _j, const FIELD _value):i(_i),j(_j),value(_value){}
  template<class ITERATOR> 
  IndexValueTriple(const ITERATOR& it):i(it.i()),j(it.j()),value(*it){}
  string str(){ostringstream result; result<<"("<<i<<","<<j<<","<<value<<")"; return result.str();}
  INDEX i; INDEX j; FIELD value; 
};

struct BlockIndexPair{
  BlockIndexPair(): block(-1), index(-1){}
  BlockIndexPair(const INDEX& _block, const INDEX& _index): block(_block), index(_index){}
  string str(){ostringstream result; result<<"("<<block<<","<<index<<")"; return result.str();}
  INDEX block; INDEX index; };

class IndexSet{
public:
  IndexSet(const int _k):k(_k){ix=new INDEX[k];}
  INDEX& operator[](const int i){return ix[i];}
  INDEX operator[](const int i) const {return ix[i];}
  void sort(){std::sort(ix,ix+k);}
  string str() const{ostringstream result; result<<"("; for(int i=0; i<k-1; i++) result<<ix[i]<<","; result<<ix[k-1]<<")"; return result.str();}
  int k;
  int* ix;
};

typedef vector<int> BlockStructure; 

struct SVpair{
  SVpair()=default;
  SVpair(const INDEX& _first, const FIELD& _second):first(_first),second(_second){}
  string str(){ostringstream result; result<<"("<<first<<","<<second<<")"; return result.str();}  
  INDEX first; 
  FIELD second;
};

struct{bool operator()(const SVpair& a, const SVpair& b){return a.first<b.first;} } SVpairComparator;

struct CSCmatrix{
  CSCmatrix()=default;
  CSCmatrix(INDEX* _ir, INDEX* _jc, FIELD* _val, const int &_nnz, const int &_nrows, const int &_ncols): ir(_ir), jc(_jc), val(_val), nnz(_nnz), nrows(_nrows), ncols(_ncols){}
  INDEX *ir, *jc;
  FIELD *val;
  int nnz, nrows, ncols;
};

class CoutLock{
public:
  CoutLock(): lock(mx){}
  lock_guard<mutex> lock;
  static mutex mx;
};


// Helper functions -------------------------------------------------------------------------------


template<typename TYPE>
inline void swapp(TYPE& x, TYPE& y) {TYPE t=x; x=y; y=t;}

template<typename TYPE>
inline void replace(TYPE*& ptr, TYPE* result) {delete ptr; ptr=result;}

template<typename TYPE>
inline void snatch(TYPE*& dest, TYPE*& source) {dest=source; source=nullptr;}

template<typename TYPE>
inline void move_over(TYPE*& ptr, TYPE*& result) {delete ptr; ptr=result; result=nullptr;}

template<typename TYPE> // this is what std::move is for 
inline TYPE&& pilfer(TYPE& x){return reinterpret_cast<TYPE&&>(x);}

template<typename TYPE>
inline TYPE movedelete(TYPE* xptr){
  TYPE x(move(*xptr)); delete xptr; return move(x);}

template<class T1, class T2>
inline bool pairLess(const  pair<T1,T2>& a, const pair<T1,T2>& b){return a.second<b.second;} 

//class pairLess{
//  bool operator()(const pair<T1,T2>& a, const pair<T1,T2>& b){return a.second<b.second;} 
//};

template<class T1, class T2>
class pairMore{
  bool operator()(const pair<T1,T2>& a, const pair<T1,T2>& b){return a.second>b.second;} 
};



// Abbreviations ----------------------------------------------------------------------------------

typedef pair<int,int> ipair;
typedef pair<double,double> dpair;

typedef MatrixX<Vectorv> MatrixXv;
typedef MatrixX<Vectorl> MatrixXl;
typedef MatrixX<Vectorh> MatrixXh;

template<class VECTOR> using Bvector=BlockedVector<VECTOR>; 
typedef BlockedVector<Cvector> BlockedCvector;
typedef BlockedVector<Vectorv> BlockedVectorv;
typedef BlockedVector<Vectorl> BlockedVectorl;
typedef BlockedVector<Vectorh> BlockedVectorh;

template<class MATRIX> using Bmatrix=BlockedMatrix<MATRIX>; 
typedef BlockedMatrix<Cmatrix> BlockedCmatrix;
typedef BlockedMatrix<MatrixX<Vectorv> > BlockedMatrixXv;
typedef BlockedMatrix<MatrixX<Vectorl> > BlockedMatrixXl;
typedef BlockedMatrix<MatrixX<Vectorh> > BlockedMatrixXh;

typedef BlockStructure Bstructure;

typedef MatrixXv MSSBLOCK;

template<class TYPE>
void debug(TYPE v){cout<<v<<endl;}

template<class TYPE1, class TYPE2>
void debug(TYPE1 v1, TYPE2 v2){cout<<v1<<v2<<endl;}

template<class TYPE1, class TYPE2, class TYPE3>
void debug(TYPE1 v1, TYPE2 v2, TYPE3 v3){cout<<v1<<v2<<v3<<endl;}

//template<class TYPE>
//void debug(string msg, TYPE v){cout<<msg<<v<<endl;}

template<class TYPE>
TYPE squared(TYPE v){return v*v;}

// Headers -----------------------------------------------------------------------------------------


#include "Rstream.hpp"
#include "Serializable.hpp"
#include "Log.hpp"



#endif


