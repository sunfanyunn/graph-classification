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


#ifndef _MatrixX
#define _MatrixX

#include <random>

#include "Matrix.hpp"
#include "Cmatrix.hpp"
#include "Cvector.hpp"
#include "Vectorv.hpp"
#include "Vectorl.hpp"
#include "Vectorh.hpp"
#include "MatrixIF.hpp"
#include "MatrixOF.hpp"

#ifdef _withEigen 
#include "EigenInterface.hpp"
#endif

extern default_random_engine randomNumberGenerator;

template<class COLUMNTYPE>
class MatrixX: public SparseMatrix, public Serializable{
public:

  MatrixX(const MatrixX<COLUMNTYPE>& x);  
  MatrixX(MatrixX<COLUMNTYPE>&& x);  
  MatrixX<COLUMNTYPE>& operator=(const MatrixX<COLUMNTYPE>& x);
  MatrixX<COLUMNTYPE>& operator=(MatrixX<COLUMNTYPE>&& x);
  MatrixX<COLUMNTYPE> copy();
  ~MatrixX(){for(int i=0; i<column.size(); i++) delete column[i];}


public: // constructors

  MatrixX(): MatrixX(0,0){};
  MatrixX(const int _nrows, const int _ncols): SparseMatrix(_nrows,_ncols), column(ncols){
    for(int i=0; i<ncols; i++) column[i]=new COLUMNTYPE(nrows);}
  //MatrixX(const int _nrows, const int _ncols, const Zero& dummy): MatrixX(_nrows,_ncols){}
  //MatrixX(const int _nrows, const int _ncols, const Uninitialized& dummy): SparseMatrix(_nrows,_ncols){}
  MatrixX(const int _nrows, const int _ncols, const Nullcols& dummy): 
    SparseMatrix(_nrows,_ncols), column(ncols,nullptr){}
  
  MatrixX(const int _nrows, const int _ncols, const Random& dummy); // depricated
  MatrixX(const int _nrows, const class Identity& dummy); // obsolete 
  MatrixX(const int _nrows, const class Random& dummy); // obsolete 


public: // named constructors

  static MatrixX Zero(const int _nrows, const int _ncols){return MatrixX(_nrows,_ncols);}
  static MatrixX Identity(const int n);
  static MatrixX Random(const int _nrows, const int _ncols, const double p=0.5);
  static MatrixX RandomSymmetric(const int n, const double p=0.5);


public: // conversions

  template<class COLUMNTYPE2>
  MatrixX(const MatrixX<COLUMNTYPE2>& x);
  template<class COLUMNTYPE2>
  MatrixX(const MatrixX<COLUMNTYPE2>& x, const int _nrows, const int _ncols);
  template<class COLUMNTYPE2>
  MatrixX(const MatrixX<COLUMNTYPE2>& x, const Transpose& dummy);
  template<class COLUMNTYPE2>
  MatrixX(const MatrixX<COLUMNTYPE2>& x, const Remap& remap1, const Remap& remap2):
    MatrixX(x,remap1,remap2,x.nrows,x.ncols){}
  MatrixX(const Cmatrix& x);

  //Cmatrix Cmatrix() const;
  //MatrixX<Vectorv> MatrixXv() const;
  //MatrixX<Vectorl> MatrixXl() const;
  //MatrixX<Vectorh> MatrixXh() const;
  


#ifdef _withEigen 
  EigenSparseMatrix convertToEigen() const;
  operator EigenSparseMatrix() const;
#endif

private: // remapping

  template<class COLUMNTYPE2>
  MatrixX(const MatrixX<COLUMNTYPE2>& x, const Remap& remap1, const Remap& remap2, const int _nrows, const int _ncols);


public: // element access and views
  
  FIELD& operator()(const int i, const int j) {assert(i<nrows); assert(j<ncols); return (*column[j])(i);} 
  FIELD operator()(const int i, const int j) const {assert(i<nrows); assert(j<ncols); return (*column[j])(i);} 
  FIELD read(const int i, const int j) const {assert(i<nrows); assert(j<ncols); return column[j]->read(i);} 

  COLUMNTYPE& vcolumn(const int j){assert(j<ncols); return *column[j];}

  void (foreach)(std::function<void(const INDEX, const INDEX, FIELD&)> lambda){
    for(int j=0; j<ncols; j++) for(auto& p:*column[j]) lambda(p.first,j,p.second);}

  void (foreach)(std::function<void(const INDEX,const INDEX,const FIELD)> lambda) const{
    for(int j=0; j<ncols; j++) for(auto& p:*column[j]) lambda(p.first,j,p.second);}

  void foreach_in_column(const int j, std::function<void(const INDEX, FIELD&)> lambda){
    for(auto& p:*column[j]) lambda(p.first,p.second);}

  void foreach_in_column(const int j, std::function<void(const INDEX, const FIELD)> lambda) const{
    for(auto& p:*column[j]) lambda(p.first,p.second);}

  bool isFilled(const int i, const int j)const {assert(j<ncols); return column[j]->isFilled(i);}

  int nFilled() const {int t=0; for(int j=0; j<ncols; j++) t+=column[j]->nFilled(); return t;}

  COLUMNTYPE diag() const { COLUMNTYPE diagonal(min(nrows,ncols)); 
    for(int j=0; j<min(nrows,ncols); j++) diagonal(j) = (*column[j])(j); return (const COLUMNTYPE) diagonal;}


public: // scalar valued operations 

  int nnz() const;
  FIELD norm2() const;
  FIELD diff2(const MatrixX<COLUMNTYPE>& X) const;
  FIELD spectralNorm() const;
  // INDEX argmax() const;


public: // vector valued operations 

  Cvector operator*(const Cvector& x) const{
    Cvector r = Cvector::Zero(nrows);
    assert(x.n==ncols);
    for(int j=0; j<ncols; j++){
      FIELD t=x.array[j];
      for(auto& p:*column[j])
        r(p.first)+=p.second*t;
    }
    return r;
  }

  Cvector dot(const COLUMNTYPE& v) const{
    assert(v.n==nrows);
    Cvector r(ncols);
    for(int j=0; j<ncols; j++) r.array[j]=column[j]->dot(v);
    return r;
  }


public: // matrix valued operations

  Cmatrix operator*(const MatrixX<COLUMNTYPE>& x) const{
    Cmatrix M(nrows,x.ncols);
    assert(ncols==nrows); 
    // if(this==&x) for(int i=0; i<ncols; i++) for(int j=0; j<=i; j++) {M(i,j)=column[i]->dot(*x.column[j]); M(j,i)=M(i,j);}
    // else for(int i=0; i<ncols; i++) for(int j=0; j<x.ncols; j++) M(i,j)=column[i]->dot(*x.column[j]);
    return M;
  }

  Cmatrix dot(const MatrixX<COLUMNTYPE>& x) const{
    Cmatrix M(ncols,x.ncols);
    assert(x.nrows==nrows); 
    if(this==&x) for(int i=0; i<ncols; i++) for(int j=0; j<=i; j++) {M(i,j)=column[i]->dot(*x.column[j]); M(j,i)=M(i,j);}
    else for(int i=0; i<ncols; i++) for(int j=0; j<x.ncols; j++) M(i,j)=column[i]->dot(*x.column[j]);
    return M;
  }

  Cmatrix rowGram() const{ // could be improved by exploiting sorted sparse vectors
    Cmatrix G=Cmatrix::Zero(nrows,nrows);
    for(int k=0; k<ncols; k++)
      for(auto& p:*column[k])
	for(auto& q:*column[k])
	  G(p.first,q.first)+=p.second*q.second;
    return G;
  }

  Cmatrix colGram() const{
    Cmatrix G(ncols,ncols);
    for(int i=0; i<ncols; i++) 
      for(int j=0; j<=i; j++){
	G(i,j)=column[i]->dot(*column[j]); 
	G(j,i)=G(i,j);
      }
    return G;
  }


public: // in place operations

  MatrixX<COLUMNTYPE>& operator+=(const MatrixX<COLUMNTYPE>& x){ assert(ncols==x.ncols);
    for(int j=0; j<ncols; j++) (*column[j])+=(*x.column[j]); return *this;}

  MatrixX<COLUMNTYPE>& operator-=(const MatrixX<COLUMNTYPE>& x){ assert(ncols==x.ncols);
    for(int j=0; j<ncols; j++) (*column[j])-=(*x.column[j]); return *this;}

  MatrixX<COLUMNTYPE>& multiplyRowsBy(const Cvector& v){
    for(int j=0; j<ncols; j++) (*column[j])*=v; return *this;}

  MatrixX<COLUMNTYPE>& multiplyColsBy(const Cvector& v){
    assert(v.n==ncols); for(int j=0; j<ncols; j++) (*column[j])*=v(j); return *this;}

  MatrixX<COLUMNTYPE>& divideRowsBy(const Cvector& v){
    for(int j=0; j<ncols; j++) (*column[j])/=v; return *this;}

  MatrixX<COLUMNTYPE>& divideColsBy(const Cvector& v){
    assert(v.n==ncols); for(int j=0; j<ncols; j++) (*column[j])*=(1.0/v(j)); return *this;}

  void transpose(); 

  MatrixX<COLUMNTYPE>& tidy(){
    for(int j=0; j<ncols; j++) column[j]->tidy(); return *this;}

  void symmetrize();

public: // I/O 

  string str(const Sparse dummy) const;
  string str(const Dense dummy) const {return Matrix::str(Dense());}
  string str() const{return str(Sparse());}

  //MatrixX(DenseMatrixFile& file);
  //MatrixX(SparseMatrixFile& file);

  MatrixX(MatrixIF& file);
  void saveTo(MatrixOF& file) const;
  void dump(FIELD* result);
  CSCmatrix cscformat();

  static string classname();
  MatrixX<COLUMNTYPE>(Bifstream& ifs);
  void serialize(Bofstream& ofs) const;
  void serialize(Rstream& rstream) const;


public:

  vector<COLUMNTYPE*> column; 

};



// ---- Copying ---------------------------------------------------------------------------------------------------


template<class COLUMNTYPE> 
MatrixX<COLUMNTYPE>::MatrixX(const MatrixX<COLUMNTYPE>& x): 
  SparseMatrix(x.nrows,x.ncols), column(x.ncols){
#ifdef _MATRIXCOPYWARNING
  cout<<"WARNING: MatrixX copied."<<endl;
#endif
  for(int j=0; j<ncols; j++) column[j]=new COLUMNTYPE(*x.column[j]);
}

template<class COLUMNTYPE> 
MatrixX<COLUMNTYPE>::MatrixX(MatrixX<COLUMNTYPE>&& x): 
  SparseMatrix(x.nrows,x.ncols), column(x.ncols){
#ifdef _MATRIXMOVEWARNING
  cout<<"WARNING: MatrixX moved."<<endl;
#endif
  for(int j=0; j<ncols; j++) {column[j]=x.column[j]; x.column[j]=nullptr;}
}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>& MatrixX<COLUMNTYPE>::operator=(const MatrixX<COLUMNTYPE>& x){
  nrows=x.nrows; ncols=x.ncols; for(auto p:column) delete p; column.resize(ncols); 
#ifdef _MATRIXCOPYWARNING
  cout<<"WARNING: MatrixX assigned."<<endl;
#endif
  for(int j=0; j<ncols; j++) column[j]=new COLUMNTYPE(*x.column[j]);
  return *this;
}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>& MatrixX<COLUMNTYPE>::operator=(MatrixX<COLUMNTYPE>&& x){
  nrows=x.nrows; ncols=x.ncols; for(auto p:column) delete p; column.resize(ncols); 
#ifdef _MATRIXMOVEWARNING
  cout<<"WARNING: MatrixX move-assigned."<<endl;
#endif
  for(int j=0; j<ncols; j++) {column[j]=x.column[j]; x.column[j]=nullptr;}
  return *this;
}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE> MatrixX<COLUMNTYPE>::copy(){
  MatrixX M(nrows,ncols);
  for(int j=0; j<ncols; j++) M.column[j]=column[j]->copy();
  return *this;
}


// ---- Constructors ---------------------------------------------------------------------------------------------


/*
template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(const int _nrows, const int _ncols, const Nullcols& dummy): 
  SparseMatrix(_nrows,_ncols), column(ncols,NULL){}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(const int _nrows, const int _ncols): 
  SparseMatrix(_nrows,_ncols), column(ncols){
  for(int i=0; i<ncols; i++) column[i]=new COLUMNTYPE(nrows);
}
*/

/*
template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(const int _nrows, const class Identity& dummy): 
  MatrixX(_nrows,_nrows){
  for(int i=0; i<ncols; i++) (*this)(i,i)=1;
}
*/

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE> MatrixX<COLUMNTYPE>::Identity(const int n){
  MatrixX M(n,n);
  for(int i=0; i<n; i++) M(i,i)=1;
  return M;
}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(const int _nrows, const int _ncols, const class Random& random):
  MatrixX(_nrows,_ncols){
  double p=random.p;
  uniform_int_distribution<int> distri(0,nrows-1);
  uniform_int_distribution<int> distrj(0,ncols-1);
  uniform_real_distribution<FIELD> distr(0,1);
  for(int i=0; i<p*nrows*ncols; i++)
    (*this)(distri(randomNumberGenerator),distrj(randomNumberGenerator))=distr(randomNumberGenerator);
  
}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE> MatrixX<COLUMNTYPE>::Random(const int _nrows, const int _ncols, const double p){
  MatrixX M(_nrows,_ncols);
  uniform_int_distribution<int> distri(0,_nrows-1);
  uniform_int_distribution<int> distrj(0,_ncols-1);
  uniform_real_distribution<FIELD> distr(0,1);
  for(int i=0; i<p*_nrows*_ncols; i++)
    M(distri(randomNumberGenerator),distrj(randomNumberGenerator))=distr(randomNumberGenerator);
  return M;
}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE> MatrixX<COLUMNTYPE>::RandomSymmetric(const int n, const double p){
  MatrixX M(n,n);
  uniform_int_distribution<int> distri(0,n-1);
  uniform_real_distribution<FIELD> distr(0,1);
  for(int u=0; u<p*n*n/2; u++){
    int i=distri(randomNumberGenerator);
    int j=distri(randomNumberGenerator);
    FIELD v=distr(randomNumberGenerator);
    M(i,j)=v;
    M(j,i)=v;
  }
  return M;
}


// ---- Conversions -----------------------------------------------------------------------------------------------



template<class COLUMNTYPE> 
template<class COLUMNTYPE2>
MatrixX<COLUMNTYPE>::MatrixX(const MatrixX<COLUMNTYPE2>& x): 
  SparseMatrix(x.nrows,x.ncols), column(x.ncols){
  for(int j=0; j<ncols; j++) column[j]=new COLUMNTYPE(*x.column[j]);
  cout<<"Warning: MatrixX convert-copied."<<endl;
}

template<class COLUMNTYPE> 
template<class COLUMNTYPE2>
MatrixX<COLUMNTYPE>::MatrixX(const MatrixX<COLUMNTYPE2>& x, const int _nrows, const int _ncols): 
  SparseMatrix(_nrows,_ncols), column(_ncols){
  for(int j=0; j<ncols; j++){
    column[j]=new COLUMNTYPE(nrows);
    if (j<x.ncols) {
    for(auto& p:*x.column[j])
      if(p.first<nrows) column[j]->insert(p.first,p.second);
  }
  }
}

template<class COLUMNTYPE> 
template<class COLUMNTYPE2>
MatrixX<COLUMNTYPE>::MatrixX(const MatrixX<COLUMNTYPE2>& x, const Transpose& dummy):
  SparseMatrix(x.ncols,x.nrows), column(x.nrows){
  for(int j=0; j<ncols; j++) column[j]=new COLUMNTYPE(ncols);
  for(int i=0; i<nrows; i++){
    x.column[i]->sort();
    for(auto& p:*x.column[i])
      if(p.second!=0) column[p.first]->append(i,p.second);
  }
}

template<class COLUMNTYPE>
template<class COLUMNTYPE2> // Warning: won't work for mapping unsorted v or h to l 
MatrixX<COLUMNTYPE>::MatrixX(const MatrixX<COLUMNTYPE2>& x, 
					 const Remap& remap1, const Remap& remap2, 
					 const int _nrows, const int _ncols):
  MatrixX(_nrows,_ncols){
  for(int j=0; j<x.ncols; j++){
    const COLUMNTYPE& xcol=*x.column[j];
    COLUMNTYPE& col=*column[remap2.forward[j]];
    for(auto& p:xcol)
      col.append(remap1.forward[p.first],p.second);
  }
}

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(const class Cmatrix& x): 
  SparseMatrix(x.nrows,x.ncols), column(x.ncols){
  for(int j=0; j<ncols; j++)
    column[j]=new COLUMNTYPE(x.vcolumn(j));
  cout<<"WARNING: dense->sparse conversion"<<endl;
}


/*
template<class COLUMNTYPE>
Cmatrix MatrixX<COLUMNTYPE>::Cmatrix() const {return Cmatrix(*this);}
template<class COLUMNTYPE>
MatrixX<Vectorv> MatrixX<COLUMNTYPE>::MatrixXv() const {return MatrixX<Vectorv>(*this);}
template<class COLUMNTYPE>
MatrixX<Vectorl> MatrixX<COLUMNTYPE>::MatrixXl() const {return MatrixX<Vectorl>(*this);}
template<class COLUMNTYPE>
MatrixX<Vectorh> MatrixX<COLUMNTYPE>::MatrixXh() const {return MatrixX<Vectorh>(*this);}
*/


// ---- Operations ----------------------------------------------------------------------------------------------


 
template<class COLUMNTYPE>
void MatrixX<COLUMNTYPE>::transpose(){ 
  vector<COLUMNTYPE*> newcolumn(nrows); 
  for(int i=0; i<nrows; i++) newcolumn[i]=new COLUMNTYPE(ncols);
  for(int i=0; i<ncols; i++){
    for(auto& p:*column[i])
      if(p.second!=0) newcolumn[p.first]->append(i,p.second);
    delete column[i];
  }
  column=newcolumn; swapp(nrows,ncols);
}

template<class COLUMNTYPE>
void MatrixX<COLUMNTYPE>::symmetrize() {
  for(int j=0;j<ncols;j++){
    for(auto& p: *column[j]) {
      if (p.first<=j){
        FIELD& val=(*column[p.first])(j);
        val += p.second;
        p.second = val;
      } else {
        FIELD& val=(*column[p.first])(j);
        if(val==0)
          val += p.second;
      }
    }
  }
}


// ---- Scalar valued operations ---------------------------------------------------------------------------------


//template<class COLUMNTYPE>
//INDEX MatrixX<COLUMNTYPE>::argmax() const{
//  cout<<"Unimplemented"<<endl;
//}


template<class COLUMNTYPE>
FIELD MatrixX<COLUMNTYPE>::norm2() const{
  FIELD t=0; for(int j=0; j<ncols; j++) t+=column[j]->norm2(); return t;}


template<class COLUMNTYPE>
FIELD MatrixX<COLUMNTYPE>::diff2(const MatrixX<COLUMNTYPE>& X) const{
  assert(X.ncols==ncols);
  FIELD t=0; for(int j=0; j<ncols; j++) t+=column[j]->diff2(*X.column[j]);
  return t;
}

template<class COLUMNTYPE>
int MatrixX<COLUMNTYPE>::nnz() const{
  int result=0;
  for(int j=0;j<ncols;j++) result+=column[j]->nnz();
  return result;
}

template<class COLUMNTYPE>
FIELD MatrixX<COLUMNTYPE>::spectralNorm() const{  
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

//template<class COLUMNTYPE>
//bool MatrixX<COLUMNTYPE>::isSparse() const {
//	return true;}


// ---- I/O --------------------------------------------------------------------------------------------------------


/*
template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(SparseMatrixFile& file): 
  MatrixX(file.nrows,file.ncols){
  IndexValueTriple p; file>>p;
  while(p.i!=-1){column[p.j]->insert(p.i,p.value); file>>p;}
  tidy();
  //for(auto it=file.begin(); it!=file.end(); ++it)
  //column[(*it).j]->insert((*it).i,(*it).value);
}


template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(DenseMatrixFile& file): 
  MatrixX(file.nrows,file.ncols){
  auto it=file.begin();
  for(int i=0; i<nrows; i++)
    for(int j=0; j<ncols; j++)
      {FIELD v; file>>v; column[j]->append(i,v);}
      //{column[j]->append(i,*it); ++it;}
}
*/


template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(MatrixIF& file): MatrixX<COLUMNTYPE>(file.nrows,file.ncols){
  file.rewind();
  IndexValueTriple t;
  file>>t;
  while(t.i>=0){
    column[t.j]->insert(t.i,t.value);
    file>>t;
  }
}


template<class COLUMNTYPE>
void MatrixX<COLUMNTYPE>::saveTo(MatrixOF& file) const{
  assert(file.nrows==nrows);
  assert(file.ncols==ncols);
  if(file.sparse){
    (foreach)([&file](const int i, const int j, const FIELD v){file<<IndexValueTriple(i,j,v);});
  }else{
    for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
	file<<read(i,j);
  }
}


template<class VECTOR>
string MatrixX<VECTOR>::classname(){return "MatrixX<"+VECTOR::classname()+">";}


template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(Bifstream& ifs):SparseMatrix(0,0){
  ifs.check(classname().c_str(),0);
  ifs.read(nrows);
  ifs.read(ncols);
  for(int j=0; j<ncols; j++)
    column.push_back(new COLUMNTYPE(ifs));
};


template<class COLUMNTYPE>
void MatrixX<COLUMNTYPE>::serialize(Bofstream& ofs) const{
  ofs.tag(classname().c_str(),0);
  ofs.write(nrows);
  ofs.write(ncols);
  for(int j=0; j<ncols; j++) 
    column[j]->serialize(ofs);
};


template<class COLUMNTYPE>
void MatrixX<COLUMNTYPE>::serialize(Rstream& rstream) const{
  rstream<<"MatrixX{"<<Rstream::endl;
  rstream.var("nrows",nrows);
  rstream.var("ncols",ncols);
  for(int j=0; j<ncols; j++){
    rstream<<"  column["<<j<<"]=";
    column[j]->serialize(rstream);
  }
  rstream<<"}"<<Rstream::endl;

};


template<class COLUMNTYPE>
string MatrixX<COLUMNTYPE>::str(const Sparse dummy) const{
  ostringstream stream;
  stream.precision(5);
  for(int j=0; j<ncols; j++) 
    for(auto& it: *column[j]) stream<<"("<<it.first<<","<<j<<") : "<<it.second<<endl;
      return stream.str();
  }

template<class COLUMNTYPE>
CSCmatrix MatrixX<COLUMNTYPE>::cscformat() {
 CSCmatrix result;
  result.nnz = nnz();
  result.ir = new INDEX[result.nnz];
  result.jc = new INDEX[ncols+1];
  result.val = new FIELD[result.nnz];
  result.nrows = nrows; result.ncols = ncols; 

  result.jc[ncols] = result.nnz;

  int counter=0;
  for(int j=0;j<ncols;j++) {
    result.jc[j]=counter;
    for (auto& it: *column[j]) {
      result.ir[counter] = it.first; 
      result.val[counter] = it.second; 
      counter++;
    }
  }
  return result;
  }

template<class COLUMNTYPE>
void MatrixX<COLUMNTYPE>::dump(FIELD* result)  {
    int counter=0;
    for (int j=0;j<ncols;j++) {
      int i=0;
      for(auto&it: *column[j]) {
        while(i!=it.first) {
          result[counter]=0; counter++;
          i++;
        }
        result[counter]=it.second; counter++;
      }
    }
  }

#ifdef _withEigen
template<class COLUMNTYPE>
EigenSparseMatrix MatrixX<COLUMNTYPE>::convertToEigen() const{
  EigenSparseMatrix M(nrows,ncols);
  vector<Eigen::Triplet<double>> triplets;
  for(int j=0;j<ncols;j++) {
    for (auto& it: *column[j]) {
      triplets.push_back(Eigen::Triplet<double>(it.first,j,it.second));
    }
  }
	M.setFromTriplets(triplets.begin(), triplets.end());
	return M;
}
#endif


#ifdef _withEigen
template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::operator EigenSparseMatrix() const{
  EigenSparseMatrix M(nrows,ncols);
  vector<Eigen::Triplet<double>> triplets;
  for(int j=0;j<ncols;j++) {
    for (auto& it: *column[j]) {
      triplets.push_back(Eigen::Triplet<double>(it.first,j,it.second));
    }
  }
  M.setFromTriplets(triplets.begin(), triplets.end());
  return M;
}
#endif




template<class COLUMNTYPE>
ostream& operator<<(ostream& stream, const MatrixX<COLUMNTYPE>& x){stream<<x.str(Sparse()); return stream;}




#endif 


/*
  //MatrixX(ASCIIfile::DenseMatrix& file);
  // MatrixX(ASCIIfile::SparseMatrix& file);

template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(ASCIIfile::DenseMatrix& file): 
  MatrixX(file.nrows,file.ncols){
  float t;
  for(int i=0; i<nrows; i++)
    for(int j=0; j<ncols; j++){
      file.ifs>>t; column[j]->append(i,t);}
}


SparseMatrixv(const TowerOfBlocks<SparseMatrixv>& T, const RemapBlock& remap);

//FIELD columnDot(const int j1, const int j2) const{ return column[j1]->dot(*column[j2]);}


DenseMatrix* dotSelf() const;
template<class COLUMNTYPE>
DenseMatrix* MatrixX<COLUMNTYPE>::dotSelf() const{
  DenseMatrix* M=newDenseMatrix(ncols,ncols);
  for(int i=0; i<ncols; i++)
    for(int j=0; j<=i; j++){
      (*M)(i,j)=column[i]->dot(*column[j]);
      (*M)(j,i)=(*M)(i,j);}
  return M;
}


template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(const string filename): 
  Matrix(0,0){
  ifstream ifs(filename.c_str());
  ifs>>nrows>>ncols;
  for(int i=0; i<ncols; i++) column.push_back(new COLUMNTYPE(nrows));
  while(ifs.good()) {
    int i; int j; float v;
    ifs>>i>>j>>v;
    (*this)(i,j)=v;
  }
  ifs.close();
}

//void addColumns(const int n);

//template<class COLUMNTYPE>
//void MatrixX<COLUMNTYPE>::addColumns(const int n){
//  for(int i=0; i<n; i++) column.push_back(new COLUMNTYPE(nrows));
//  ncols+=n;
//}


template<class COLUMNTYPE>
MatrixX<COLUMNTYPE>::MatrixX(ASCIIfile::SparseMatrix& file): 
  MatrixX(file.nrows,file.ncols){
  int i; int j; float v;
  while(file.ifs.good()) {
    file.ifs>>i>>j>>v;
    column[j]->insert(i,v);
  }
}

*/

//  MatrixX<COLUMNTYPE>* newof(){return new MatrixX<COLUMNTYPE>(*this);}

//template<class COLUMNTYPE>
//class MatrixXiterator;

/*
public: // iterator 

  typedef MatrixXiterator<COLUMNTYPE> iterator;

  iterator begin() const{
    assert(column.size()>0); 
    for(int j=0; j<ncols; j++){
      assert(column[j]!=NULL);
      if(column[j]->size()>0) return iterator(*this,column[j]->begin()->first,j,column[j]->begin());}
    return iterator(*this);
  }

  void increment(iterator& itr) const{
    itr.it++; 
    while(itr.it==column[itr._j]->end()){
      if(++itr._j>=ncols){itr.endflag=true; return;}
      itr.it=column[itr._j]->begin();
    }
    itr._i=itr.it->first; 
  }

  iterator end() const{ return iterator(*this);}
*/

/*
template<class COLUMNTYPE>
class MatrixXiterator{
  friend class MatrixX<COLUMNTYPE>;
public:
  MatrixXiterator(const MatrixX<COLUMNTYPE>& _owner):owner(_owner), endflag(true){};
  MatrixXiterator(const MatrixX<COLUMNTYPE>& _owner, 
			const int __i, const int __j, const typename COLUMNTYPE::iterator& _it): 
    owner(_owner),_i(__i),_j(__j),it(_it){}
  FIELD& operator*(){return it->second;}
  FIELD operator*()const{return it->second;}
  FIELD* operator->(){return &it->second;}
  int i()const{return _i;}
  int j()const{return _j;}
  MatrixXiterator& operator++(){owner.increment(*this); return *this;}
  void operator++(int unused){owner.increment(*this);}
  bool operator!=(const MatrixXiterator& x){return endflag!=x.endflag;}
private:
  int _i;
  int _j;
  typename COLUMNTYPE::iterator it;
  const MatrixX<COLUMNTYPE>& owner;
  bool endflag=false;
};
*/


