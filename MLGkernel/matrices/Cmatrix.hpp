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


#ifndef _Cmatrix
#define _Cmatrix

#include "Matrix.hpp"
#include "Cvector.hpp"

#ifdef _withEigen
class EigenVectorXdAdaptor;
class EigenMatrixXdAdaptor;
#endif


class Cmatrix: public DenseMatrix, public Serializable{
public:

  Cmatrix(const Cmatrix& x);
  Cmatrix(Cmatrix&& x);
  Cmatrix(AtomicCmatrix&& x);
  Cmatrix& operator=(const Cmatrix& x);
  Cmatrix& operator=(Cmatrix&& x);
  Cmatrix copy() const;
  ~Cmatrix(){delete[] array;}


public: 

  Cmatrix(): Cmatrix(0,0){array=nullptr;}
  Cmatrix(const int _nrows, const int _ncols): 
    DenseMatrix(_nrows,_ncols) {array=new FIELD[nrows*ncols];}
  Cmatrix(const int _nrows, const int _ncols, const Zero dummy): // obsolete     
    Cmatrix(_nrows,_ncols) {for(int i=0; i<nrows*ncols; i++) array[i]=0;}


public: // named constructors 

  static Cmatrix Zero(const int _nrows, const int _ncols){
    Cmatrix M(_nrows,_ncols); for(int i=0; i<_nrows*_ncols; i++) M.array[i]=0; return M;}

  static Cmatrix Identity(const int n){
    Cmatrix M=Zero(n,n); for(int i=0; i<n; i++) M.array[i*n+i]=1; return M;}

  static Cmatrix Kronecker(const Cmatrix& seed, const int nlevels);
  static Cmatrix Kronecker(const Cmatrix& X1, const Cmatrix& X2);

  static Cmatrix Random(const int _nrows, const int _ncols);
  static Cmatrix RandomSymmetric(const int _nrows);
  static Cmatrix GaussianSymmetric(const int _nrows);
  static Cmatrix RandomPosDef(const int n);


public: // converters

  Cmatrix(const Cmatrix& x, const Remap& remap1, const Remap& remap2); 
  Cmatrix(const Cmatrix& x, const Transpose& dummy);

  template<class COLUMNTYPE> 
  Cmatrix(const MatrixX<COLUMNTYPE>& x); 

  //Cmatrix(const int _nrows, const int _ncols, const FIELD* _array);

#ifdef _withEigen
  Cmatrix(const EigenMatrixXdAdaptor& M);
  Cmatrix(const EigenMatrixXdAdaptor& M, const int _nrows, const int _ncols=0);
#endif  

  template<class TYPE> TYPE convert() const; 


public: // Comparisons

  bool operator==(const Cmatrix& X) const{ 
    if(X.nrows!=nrows) return false;
    if(X.ncols!=ncols) return false;
    for(int i=0; i<nrows*ncols; i++) if(array[i]!=X.array[i]) return false;
    return true;
  }

  bool operator!=(const Cmatrix& X) const {return !((*this)==X);}


public: // element access and views
  
  FIELD& operator()(const int i, const int j){
    assert(i<nrows); assert(j<ncols); return array[j*nrows+i];} 
  FIELD operator()(const int i, const int j) const{ 
    assert(i<nrows); assert(j<ncols); return array[j*nrows+i];} 
  FIELD read(const int i, const int j) const{ 
    assert(i<nrows); assert(j<ncols); return array[j*nrows+i];} 

  void (foreach)(std::function<void(const INDEX, const INDEX,FIELD&)> lambda){
    for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) lambda(i,j,array[j*nrows+i]);}

  void (foreach)(std::function<void(const INDEX, const INDEX, const FIELD)> lambda) const{
    for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) lambda(i,j,array[j*nrows+i]);}

  void foreach_in_column(const int j, std::function<void(const INDEX, FIELD&)> lambda){
    for(int i=0; i<nrows; i++) lambda(i,array[j*nrows+i]);}

  void foreach_in_column(const int j, std::function<void(const INDEX, const FIELD)> lambda) const{
    for(int i=0; i<nrows; i++) lambda(i,array[j*nrows+i]);}

  Cvector::Virtual vcolumn(const int j) const {
    return Cvector::Virtual(nrows,&array[j*nrows]);}

  Cvector row(const int i) const{
    Cvector v(ncols); for(int j=0; j<ncols; j++) v(j)=array[j*nrows+i]; return v;}

  Cvector column(const int j) const{
    Cvector v(nrows); for(int i=0; i<nrows; i++) v(i)=array[j*nrows+i]; return v;}

  Cmatrix submatrix(const int _nrows, const int _ncols, const int ioffset=0, const int joffset=0) const;
  Cmatrix submatrix(const IndexSet& I1, const IndexSet& I2) const{
    const int k1=I1.k;
    const int k2=I2.k;
    Cmatrix R(k1,k2);
    for(int i=0; i<k1; i++)
      for(int j=0; j<k2; j++)
	R(i,j)=array[I2[j]*nrows+I1[i]];
    return R;
  }

  Cvector diag() const {Cvector diagonal(min(nrows,ncols)); 
    for(int j=0; j<min(nrows,ncols); j++) diagonal(j) = array[j*(nrows+1)]; return diagonal;}


public: // scalar-valued operations

  int nnz() const {int t=0; for(int i=0; i<nrows*ncols; i++) if(array[i]!=0) t++; return t;}

  FIELD norm2() const{
    FIELD t=0; for(int i=0; i<nrows*ncols; i++) t+=array[i]*array[i]; return t;}

  FIELD diff2(const Cmatrix& X) const{
    assert(X.nrows==nrows); assert(X.ncols==ncols); FIELD t=0;
    for(int i=0; i<nrows*ncols; i++) t+=(array[i]-X.array[i])*(array[i]-X.array[i]);
    return t;}

    
  FIELD spectralNorm() const{  
    FIELD norm=0;
    Cvector v = Cvector::Random(ncols);
    FIELD oldnorm=0;
    
    for(int i=0; i<200; i++){
      cout<<"test"<<endl;
      v=(*this)*v; //-M*v;
      norm=v.norm2();
      if(fabs((norm-oldnorm)/norm)<0.001) break;
      oldnorm=norm;
      v*=1.0/sqrt(norm);
    }
    return norm;
    
  }
    
  FIELD columnsum(const int j) const{
    FIELD t=0; for(int i=0; i<nrows; i++) t+=array[j*nrows+i]; return t;}


public: // vector valued operations 

  Cvector operator*(const Cvector& x) const{
    Cvector y=Cvector::Zero(nrows);
    for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
	       y.array[i]+=array[j*nrows+i]*x.array[j];
    return y;
  }
  
  Cvector dot(const Cvector& v) const{
    assert(v.n==nrows);
    Cvector r(ncols);
    for(int j=0; j<ncols; j++){
      FIELD t=0; for(int k=0; k<nrows; k++) t+=array[j*nrows+k]*v.array[k];
      r.array[j]=t;}
    return r;
  }

    
public: // matrix valued operations

  Cmatrix operator*(const FIELD x) const{
    Cmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]*x; return M;}

  Cmatrix operator*(const Cmatrix& X) const{
    Cmatrix M=Zero(nrows,X.ncols);
    assert(ncols==X.nrows);
    for(int i=0; i<nrows; i++)
      for(int j=0; j<X.ncols; j++)
	for(int k=0; k<ncols; k++)
	  M.array[j*nrows+i]+=array[k*nrows+i]*X.array[j*X.nrows+k];
    return M;
  }
  
  Cmatrix operator+(const Matrix& x) const{
    Cmatrix R(nrows,ncols);
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
	R(i,j)=(*this)(i,j)+x(i,j);
    return R;
  }

 Cmatrix operator-(const Matrix& x) const{
    Cmatrix R(nrows,ncols);
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
	R(i,j)=(*this)(i,j)-x(i,j);
    return R;
  }

  Cmatrix dot(const Cmatrix& x) const{
    Cmatrix M(ncols,x.ncols);
    if(this==&x) for(int i=0; i<ncols; i++) for(int j=0; j<=i; j++){
      FIELD t=0; for(int k=0; k<nrows; k++) t+=array[i*nrows+k]*x.array[j*x.nrows+k]; M(i,j)=t; M(j,i)=t;}
    else for(int i=0; i<ncols; i++) for(int j=0; j<x.ncols; j++){
      FIELD t=0; for(int k=0; k<nrows; k++) t+=array[i*nrows+k]*x.array[j*x.nrows+k]; M(i,j)=t;}
    return M;
  }
  
  Cmatrix rowGram() const{
    Cmatrix G(nrows,nrows);
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<=i; j++){
	FIELD t=0; 
	for(int k=0; k<ncols; k++) 
	  t+=array[k*nrows+i]*array[k*nrows+j]; 
	G(i,j)=t; G(j,i)=t;
      }
    return G;
  }

  Cmatrix colGram() const{
    Cmatrix G(ncols,ncols);
    for(int i=0; i<ncols; i++) 
      for(int j=0; j<=i; j++){
	FIELD t=0; 
	for(int k=0; k<nrows; k++)
	  t+=array[i*nrows+k]*array[j*nrows+k]; 
	G(i,j)=t; G(j,i)=t;
      }
    return G;
  }

  //Cmatrix colGram(const Cmatrix& dummy) const {return colGram();}


public: // in-place operations 


  Cmatrix& operator*=(const FIELD& x){
    for(int i=0; i<nrows*ncols; i++) array[i]*=x; return *this;}


  Cmatrix& operator+=(const Cmatrix& x){
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    for(int i=0; i<nrows*ncols; i++) array[i]+=x.array[i]; return *this;}


  Cmatrix& operator-=(const Cmatrix& x){
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    for(int i=0; i<nrows*ncols; i++) array[i]-=x.array[i]; return *this;}


  Cmatrix& multiplyRowsBy(const Cvector& v){
    assert(v.n==nrows);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	array[j*nrows+i]=v.array[i]*array[j*nrows+i];
    return *this;}


  Cmatrix& multiplyColsBy(const Cvector& v){
    assert(v.n==ncols);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	array[j*nrows+i]=v.array[j]*array[j*nrows+i];
    return *this;}

  Cmatrix& divideRowsBy(const Cvector& v){
    assert(v.n==nrows);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	array[j*nrows+i]=array[j*nrows+i]/v.array[i];
    return *this;}

  Cmatrix& divideColsBy(const Cvector& v){
    assert(v.n==ncols);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	array[j*nrows+i]=array[j*nrows+i]/v.array[j];
    return *this;}

  void symmetrize();


public: // solvers

  FIELD determinant() const;
  Cmatrix exp() const;
  Cmatrix inverse() const;
  Cmatrix inverseSqrt() const;
  Cvector solve(const Cvector& b) const;
  pair<Cmatrix*,Cvector*> symmetricEigensolver() const;
  triple<Cmatrix,Cmatrix,Cvector> svd() const;
  Cmatrix exp(const FIELD t) const;

  Cmatrix eigenvectors(const int k=0);
  Cvector eigenvalues(const int k=0) const;

public: // I/O

  string str(const Dense dummy) const; 
  string str(const Sparse dummy) const;
  string str() const; 

  //Cmatrix(DenseMatrixFile& file);
  //Cmatrix(SparseMatrixFile& file);

  Cmatrix(MatrixIF& file);
  void saveTo(MatrixOF& file) const;
  void dump(FIELD* result);
  CSCmatrix cscformat();

  static string classname();
  Cmatrix(Bifstream& ifs);
  void serialize(Bofstream& ofs) const;
  void serialize(Rstream& rstream) const;


public: // iterators

  //CmatrixIterator begin() const; 
  //void increment(CmatrixIterator& itr) const; 
  //CmatrixIterator end() const; 


public: // depricated 

  // Cmatrix(const int _nrows, const Identity dummy):                   
  //   Cmatrix(_nrows,_nrows,Zero()) {for(int i=0; i<nrows; i++) array[i*nrows+i]=1;}
  // Cmatrix(const int _nrows, const int _ncols, const Random dummy);
  // Cmatrix(const int _nrows, const Random dummy);
  // Cmatrix* newof(){return new Cmatrix(*this);}

  //  Cvector::Virtual column(const int j) const {
  //  return Cvector::Virtual(nrows,&array[j*nrows]);}


  
public:

  FIELD* array;

};


ostream& operator<<(ostream& stream, const Cmatrix& v);



// ---- Template methods -----------------------------------------------------------------------------------------


template<class COLUMNTYPE>
Cmatrix::Cmatrix(const MatrixX<COLUMNTYPE>& x): 
  Cmatrix(x.nrows,x.ncols){
  for(int i=0; i<nrows*ncols; i++) array[i]=0;
  cout<<"Warning!!!!"<<endl;
  for(int j=0; j<ncols; j++) for(auto& p:*x.column[j]) array[j*nrows+p.first]=p.second;}



// ---- Hash function ---------------------------------------------------------------------------------------------


namespace std{
  template<>
  class hash<Cmatrix>{
  public:
    size_t operator()(const Cmatrix &M ) const{
      size_t h=0; for(int i=0; i<M.nrows*M.ncols; i++) h=(h<<1)^hash<FIELD>()(M.array[i]); return h;}
  };
};



#endif 
