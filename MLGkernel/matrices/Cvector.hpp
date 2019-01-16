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


#ifndef _Cvector
#define _Cvector

#include "DenseVector.hpp"
#include "Remap.hpp"

#ifdef _withEigen
class EigenVectorXdAdaptor;
#endif


class Cvector: public DenseVector, public Serializable{
public: 

  class Virtual;

  Cvector(const Cvector& x): Cvector(x.n){
    #ifdef _MATRIXCOPYWARNING
      cout<<"WARNING: Cvector copied."<<endl;
    #endif
    for(int i=0; i<n; i++) array[i]=x.array[i]; // {std::copy(x.array,x.array+n,array);} 
  }

  Cvector copy() const{ // silent copies 
    Cvector v(n); for(int i=0; i<n; i++) v.array[i]=array[i]; return v;}

  Cvector(Cvector&& x): DenseVector(x.n){
    #ifdef _MATRIXMOVEWARNING
      cout<<"WARNING: Cvector moved."<<endl;
    #endif
    array=x.array; x.array=nullptr;
  }

  Cvector& operator=(const Cvector& x){
    #ifdef _MATRIXCOPYWARNING
      cout<<"WARNING: Cvector assigned."<<endl;
    #endif
      n=x.n; delete[] array; array=new FIELD[n]; //std::copy(x.array,x.array+n,array); 
    for(int i=0; i<n; i++) array[i]=x.array[i]; return *this;
  }

  Cvector& operator=(Cvector&& x){
    #ifdef _MATRIXMOVEWARNING
      cout<<"WARNING: Cvector move-assigned."<<endl;
    #endif
      n=x.n; delete[] array; array=x.array; x.array=nullptr; return *this;
  }

  ~Cvector(){delete[] array;}


public: // constructors

  Cvector(): DenseVector(0){array=nullptr;}
  Cvector(const int _n): DenseVector(_n) {array=new FIELD[n];}
  // Cvector(const int _n, const Random& dummy); // deprecated
  Cvector(const int _n, const FIELD* _array);  
  Cvector(const initializer_list<FIELD> list);


public: // named constructors

  static Cvector Zero(const int n){Cvector v(n); for(int i=0; i<n; i++) v.array[i]=0; return v;}
  static Cvector Filled(const int n, const FIELD t){Cvector v(n); for(int i=0; i<n; i++) v.array[i]=t; return v;}
  static Cvector Random(const int n);
  static Cvector Remap(const Cvector& x, const class Remap& remap) {return Cvector(x,remap,false);}
  static Cvector InverseRemap(const Cvector& x, const class Remap& remap) {return Cvector(x,remap,true);}


public: // converters 

  Cvector(const Cvector& x, const class Remap& remap, const bool inverse=false): Cvector(x.n){
    if(!inverse) for(int i=0; i<n; i++) array[i]=x.array[remap.forward[i]];
    else for(int i=0; i<n; i++) array[i]=x.array[remap.backward[i]];}
  //Cvector(const Vectorv& x);
  //Cvector(const Vectorl& x);
  //Cvector(const Vectorh& x);

#ifdef _withEigen
  Cvector(const EigenVectorXdAdaptor& x);
#endif


  template<class TYPE> TYPE convert() const;

  static Cvector merge(const Cvector& x, const Cvector& y);


public: // Comparisons

  bool operator==(const Cvector& x) const{
    if(x.n!=n) return false;
    for(int i=0; i<n; i++) if(array[i]!=x.array[i]) return false;
    return true;
  }

  bool operator!=(const Cvector& x) const {return !((*this)==x);}


public: // element access 

  FIELD& operator()(const int i){return array[i];}
  FIELD operator()(const int i) const {return array[i];}

  void (foreach)(std::function<void(const INDEX, FIELD&)> lambda) {for(int i=0; i<n; i++) lambda(i,array[i]);}
  void (foreach)(std::function<void(const INDEX, const FIELD)> lambda) const {for(int i=0; i<n; i++) lambda(i,array[i]);}

  void for_each(std::function<void(const INDEX, FIELD&)> lambda) {for(int i=0; i<n; i++) lambda(i,array[i]);}
  void for_each(std::function<void(const INDEX, const FIELD)> lambda) const {for(int i=0; i<n; i++) lambda(i,array[i]);}

  int size() const {return n;}

  Cvector::Virtual vsubvector(const int beg, const int end);

public: // scalar-valued operations 

  //bool operator==(const Cvector& x) {assert(x.n==n);
  //  for(int i=0; i<n; i++) if (array[i]!=x.array[i]) return false; return true;}

  int nnz() const{int t=0; for(int i=0; i<n; i++) if (array[i]!=0) t++; return t;}

  FIELD sum() const{FIELD t=0; for(int i=0; i<n; i++) t+=array[i]; return t;}

  FIELD max() const{
    FIELD t=array[0]; for(int i=1; i<n; i++) if(array[i]>t) t=array[i]; return t;}
  FIELD max_abs() const{
    FIELD t=fabs(array[0]); for(int i=1; i<n; i++) if(fabs(array[i])>t) t=array[i]; return t;}
  INDEX argmax() const{ 
    if(n==0) return 0; int best=0; FIELD max=array[best];  
    for(int i=0; i<n; i++) if(array[i]>max) {best=i; max=array[i];}
    return best;}
  INDEX argmax_abs() const{
    if(n==0) return 0; int best=0; FIELD max=fabs(array[best]);  
    for(int i=0; i<n; i++) if(fabs(array[i])>max) {best=i; max=fabs(array[i]);}
    return best;}

  FIELD min() const{
    FIELD t=array[0]; for(int i=1; i<n; i++) if(array[i]<t) t=array[i]; return t;}
  FIELD min_abs() const{
    FIELD t=fabs(array[0]); for(int i=1; i<n; i++) if(fabs(array[i])<t) t=array[i]; return t;}
  INDEX argmin() const{ 
    if(n==0) return 0; int best=0; FIELD min=array[best];  
    for(int i=0; i<n; i++) if(array[i]<min) {best=i; min=array[i];}
    return best;}
  INDEX argmin_abs() const{ 
    if(n==0) return 0; int best=0; FIELD min=fabs(array[best]);  
    for(int i=0; i<n; i++) if(fabs(array[i])<min) {best=i; min=array[i];}
    return best;}



  FIELD norm2() const {FIELD t=0; for(int i=0; i<n; i++) t+=array[i]*array[i]; return t;}
  
  FIELD norm() const { return sqrt(norm2()); }

  FIELD diff2(const Cvector& x) const{
    assert(x.n==n); FIELD t=0; for(int i=0; i<n; i++) t+=(array[i]-x.array[i])*(array[i]-x.array[i]); return t;}

  FIELD dot(const Cvector& x) const{assert(x.n==n);
    FIELD t=0; for(int i=0; i<n; i++) t+=array[i]*x.array[i]; return t;}


public: // in-place operations 

  Cvector& operator+=(const FIELD& x) {for(int i=0; i<n; i++) array[i]+=x; return *this;}
  Cvector& operator-=(const FIELD& x) {for(int i=0; i<n; i++) array[i]-=x; return *this;}
  Cvector& operator*=(const FIELD& x) {for(int i=0; i<n; i++) array[i]*=x; return *this;}
  Cvector& operator/=(const FIELD& x) {for(int i=0; i<n; i++) array[i]/=x; return *this;}

  Cvector& operator+=(const Cvector& x) {assert(x.n==n); for(int i=0; i<n; i++) array[i]+=x.array[i]; return *this;}
  Cvector& operator-=(const Cvector& x) {assert(x.n==n); for(int i=0; i<n; i++) array[i]-=x.array[i]; return *this;}
  Cvector& operator*=(const Cvector& x) {assert(x.n==n); for(int i=0; i<n; i++) array[i]*=x.array[i]; return *this;}
  Cvector& operator/=(const Cvector& x) {assert(x.n==n); for(int i=0; i<n; i++) array[i]/=x.array[i]; return *this;}

  Cvector& add(const Cvector& x) {assert(n==x.n); for(int i=0; i<n; i++) array[i]+=x.array[i]; return *this;}
  Cvector& add(const Cvector& x, const FIELD c) {assert(n==x.n); for(int i=0; i<n; i++) array[i]+=c*x.array[i]; return *this;}


public: // vector-valued operations 

  Cvector operator*(const FIELD& x) const{
    Cvector v(n); for(int i=0; i<n; i++) v.array[i]=array[i]*x; return v;}

  Cvector operator/(const FIELD& x) const{
    Cvector v(n); for(int i=0; i<n; i++) v.array[i]=array[i]/x; return v;}

  Cvector operator-(const Cvector& x) const{
    assert(x.n==n); Cvector v(n); for(int i=0; i<n; i++) v.array[i]=array[i]-x.array[i]; return v;}

public: // I/O

  Cvector(Bifstream& ifs);
  void serialize(Bofstream& ofs) const;
  void serialize(Rstream& rstream) const;
  static string classname();


public:  

  FIELD* array;

};



class Cvector::Virtual : public Cvector{
public:
  Virtual(const int _n, FIELD* _array): Cvector(){n=_n; array=_array;}
  Virtual(const Cvector::Virtual& x): Virtual(x.n,x.array){} 
  ~Virtual(){array=NULL;}
};



ostream& operator<<(ostream& stream, const Cvector& v);


// ---- Hash function ---------------------------------------------------------------------------------------------


namespace std{
  template<>
  class hash<Cvector>{
  public:
    size_t operator()(const Cvector& v) const{
      size_t h=0; for(int i=0; i<v.n; i++) h=(h<<1)^hash<FIELD>()(v.array[i]); return h;}
  };
};



#endif
