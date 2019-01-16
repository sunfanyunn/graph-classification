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


#include "Cmatrix.hpp"
#include "AtomicCmatrix.hpp"
#include "MatrixX.hpp"
//#include "DenseMatrixFile.hpp"
//#include "SparseMatrixFile.hpp"
#include "MatrixIF.hpp"
#include "MatrixOF.hpp"

#ifdef _withEigen
#include "EigenInterface.hpp"
#endif

#ifdef _withLapack
#include "LapackInterface.hpp"
#endif

extern default_random_engine randomNumberGenerator;


// ---- Copying ---------------------------------------------------------------------------------------------------


Cmatrix::Cmatrix(const Cmatrix& x):
  Cmatrix(x.nrows,x.ncols){
  //std::copy(x.array,x.array+nrows*ncols,array);
  for(int i=0; i<nrows*ncols; i++) array[i]=x.array[i];
  #ifdef _MATRIXCOPYWARNING
    cout<<"WARNING: Cmatrix copied."<<endl;
  #endif
}

Cmatrix::Cmatrix(Cmatrix&& x):
  DenseMatrix(x.nrows,x.ncols){
  array=x.array; x.array=nullptr; x.nrows=0; x.ncols=0;
  #ifdef _MATRIXMOVEWARNING
    cout<<"WARNING: Cmatrix moved."<<endl;
  #endif
}

Cmatrix::Cmatrix(AtomicCmatrix&& x):
  DenseMatrix(x.nrows,x.ncols){
  array=x.array; x.array=nullptr; x.nrows=0; x.ncols=0;
  #ifdef _MATRIXMOVEWARNING
    cout<<"WARNING: Cmatrix moved."<<endl;
  #endif
}

Cmatrix& Cmatrix::operator=(const Cmatrix& x){
  delete[] array; nrows=x.nrows; ncols=x.ncols; array=new FIELD[nrows*ncols];
  //std::copy(x.array,x.array+nrows*ncols,array);
  for(int i=0; i<nrows*ncols; i++) array[i]=x.array[i];
  #ifdef _MATRIXCOPYWARNING
    cout<<"WARNING: Cmatrix assigned."<<endl;
  #endif
  return *this;
}

Cmatrix& Cmatrix::operator=(Cmatrix&& x){
  delete[] array;  nrows=x.nrows; ncols=x.ncols; 
  array=x.array; x.array=nullptr; x.nrows=0; x.ncols=0;
  #ifdef _MATRIXMOVEWARNING
    cout<<"WARNING: Cmatrix move-assigned."<<endl;
  #endif
  return *this;
}

Cmatrix Cmatrix::copy() const{
  Cmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]; return M;}



// ---- Constructors ----------------------------------------------------------------------------------------------

/*
Cmatrix::Cmatrix(const int _nrows, const int _ncols, const Random dummy): 
  DenseMatrix(_nrows,_ncols){
  array=new FIELD[nrows*ncols]; 
  uniform_real_distribution<FIELD> distr;
  for(int i=0; i<nrows*ncols; i++) array[i]=distr(randomNumberGenerator);
}


Cmatrix::Cmatrix(const int _nrows, const Random dummy): 
  DenseMatrix(_nrows,_nrows){
  array=new FIELD[nrows*ncols]; 
  uniform_real_distribution<FIELD> distr;
  for(int i=0; i<nrows; i++) 
    for(int j=0; j<=i; j++){
      FIELD t=distr(randomNumberGenerator); 
      array[j*nrows+i]=t; 
      array[i*nrows+j]=t;
    }
}
*/

// ---- Conversions -----------------------------------------------------------------------------------------------


Cmatrix::Cmatrix(const Cmatrix& x, const Remap& remap1, const Remap& remap2):
  Cmatrix(x.nrows,x.ncols){
  for(int i=0; i<nrows; i++) 
    for(int j=0; j<ncols; j++) 
      array[j*nrows+i]=x(remap1.forward[i],remap2.forward[j]);
}

Cmatrix::Cmatrix(const Cmatrix& x, const Transpose& dummy):
  Cmatrix(x.ncols,x.nrows){
  for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) array[j*nrows+i]=x.array[i*ncols+j];
}

/*
Cmatrix::Cmatrix(const int _nrows, const int _ncols, const FIELD* _array): Cmatrix(_nrows,_ncols){
int counter=0;
for(int j=0;j<ncols;j++)
for(int i=0;i<nrows;i++){
counter++;
(*this)(i,j) = _array[counter];
}
}
*/


#ifdef _withEigen
Cmatrix::Cmatrix(const EigenMatrixXdAdaptor& X):Cmatrix(X.rows(),X.cols()){
  for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) (*this)(i,j)=X(i,j);
}
#endif


#ifdef _withEigen
Cmatrix::Cmatrix(const EigenMatrixXdAdaptor& X, const int _nrows, const int _ncols):
  Cmatrix( (_nrows==0 ? X.rows() : _nrows), (_ncols==0 ? X.cols() : _ncols) ){
  for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) (*this)(i,j)=X(i,j);
}
#endif


#ifdef _withEigen
template<>
Eigen::MatrixXd Cmatrix::convert() const{
  Eigen::MatrixXd M(nrows,ncols);
  for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) M(i,j)=array[j*nrows+i];
  return M;
}
#endif

//MatrixXv Cmatrix::MatrixXv() const{return ::MatrixXv(*this);}
//MatrixXl Cmatrix::MatrixXl() const{return ::MatrixXl(*this);}
//MatrixXh Cmatrix::MatrixXh() const{return ::MatrixXh(*this);}




// ---- Specialized Constructors ----------------------------------------------------------------------------------



Cmatrix Cmatrix::Kronecker(const Cmatrix& seed, const int k){
  assert(k>=1);
  if(k==1) return seed;
  Cmatrix Msub=Kronecker(seed,k-1);
  return Kronecker(seed,Msub);
}


Cmatrix Cmatrix::Kronecker(const Cmatrix& X1, const Cmatrix& X2){
  int _nrows=X2.nrows*X1.nrows;
  int _ncols=X2.ncols*X1.ncols;
  Cmatrix M(_nrows,_ncols);
  for(int bJ=0; bJ<X1.ncols; bJ++)
    for(int bI=0; bI<X1.nrows; bI++){
      FIELD s=X1(bI,bJ);
      for(int j=0; j<X2.ncols; j++)
	for(int i=0; i<X2.nrows; i++)
	  M(bI*X2.nrows+i,bJ*X2.ncols+j)=s*X2(i,j);
    }
  return M;
}


Cmatrix Cmatrix::Random(const int _nrows, const int _ncols){
  Cmatrix M(_nrows,_ncols);
  uniform_real_distribution<FIELD> distr;
  for(int i=0; i<_nrows*_ncols; i++) M.array[i]=distr(randomNumberGenerator);
  return M;
}


Cmatrix Cmatrix::RandomSymmetric(const int n){
  Cmatrix M(n,n);
  uniform_real_distribution<FIELD> distr;
  for(int i=0; i<n; i++) 
    for(int j=0; j<=i; j++){
      FIELD t=distr(randomNumberGenerator); 
      M.array[j*n+i]=t; 
      M.array[i*n+j]=t;
    }
  return M;
}


Cmatrix Cmatrix::GaussianSymmetric(const int n){
  Cmatrix M(n,n);
  normal_distribution<FIELD> distr;
  for(int i=0; i<n; i++) 
    for(int j=0; j<=i; j++){
      FIELD t=distr(randomNumberGenerator); 
      M.array[j*n+i]=t; 
      M.array[i*n+j]=t;
    }
  return M;
}


Cmatrix Cmatrix::RandomPosDef(const int n){
  Cmatrix A=Random(n,n);
  return A.dot(A);
}


// ---- Submatrices -----------------------------------------------------------------------------------------------



Cmatrix Cmatrix::submatrix(const int _nrows, const int _ncols, const int ioffset, const int joffset) const{
  Cmatrix R(_nrows,_ncols);
  assert(ioffset+_nrows<=nrows);
  assert(joffset+_ncols<=ncols);
  for(int j=0; j<_ncols; j++)
    for(int i=0; i<_nrows; i++)
      R.array[j*_nrows+i]=array[(j+joffset)*nrows+i+ioffset];
  return R;
}



// ---- Binary Operations -----------------------------------------------------------------------------------------



// ---- Solvers ----------------------------------------------------------------------------------------------

FIELD Cmatrix::determinant() const{
#ifdef _withEigen
	Eigen::MatrixXd A = convert<Eigen::MatrixXd>();
	return A.determinant();
#else 
	{CoutLock lock; cout<<"Error: Cmatrix::determinant() called, but pMMF was compiled with no linear algebra library."<<endl; return 0;}
#endif
}

Cmatrix Cmatrix::exp() const{
  //#ifdef _withEigen
  //Eigen::MatrixXd A=convert<Eigen::MatrixXd>();
  //return EigenMatrixXdAdaptor(A.exp());
  //#endif
  cout<<"Exponentiation disabled"<<endl;
  return copy();
}

Cmatrix Cmatrix::inverse() const{
#ifdef _withEigen
  Eigen::MatrixXd A=convert<Eigen::MatrixXd>();
  return EigenMatrixXdAdaptor(A.inverse());
#elif _withLapack
	int IPIV[nrows];
	Cmatrix matrix_inverse(*this);
	// LU factorization
	int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR,nrows,ncols,matrix_inverse.array,nrows,IPIV);
	if(info) {CoutLock lock; cout<<"Error: DGETRF failed!"<<endl; exit(-1);}
	// Inverse using the computed LU factorization
	info = LAPACKE_dgetri(LAPACK_COL_MAJOR,nrows,matrix_inverse.array,nrows,IPIV);
	if(info) {CoutLock lock; cout<<"Error: DGETRI failed!"<<endl; exit(-1);}
	return matrix_inverse;
#else
  {CoutLock lock; cout<<"Error: Cmatrix::inverse() called, but pMMF was compiled with no linear algebra library."<<endl;}
  exit(-1);
#endif
}


Cmatrix Cmatrix::inverseSqrt() const{
#ifdef _withEigen
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(convert<Eigen::MatrixXd>());
  return EigenMatrixXdAdaptor(solver.operatorInverseSqrt());
#elif _withLapack
	Cmatrix matrix_inversesqrt(nrows,ncols);
	int IPIV[nrows];
	double W[nrows], Zprime[nrows*ncols]; 
	Cmatrix Z(*this);
	// eigendecomposition Z*W*(Z^T)
	int info = LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','L',nrows,Z.array,nrows,W);
	if(info) {CoutLock lock; cout<<"Error: DSYEV failed!"<<endl; exit(-1);}
	// Zprime = Z*inverseSqrt(diag(W))
	double temp=0; int counter=0;
	for(int j=0;j<ncols;j++){
		temp = 1.0/sqrt(W[j]);
		for(int i=0;i<nrows;i++) {
			Zprime[counter] = Z.array[counter]*temp;
			counter++;
		}
	}
	// inverse sqrt = (Zprime)*(Z^T)
	cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,nrows,nrows,nrows,1,Zprime,nrows,Z.array,nrows,0,matrix_inversesqrt.array,nrows);
	return matrix_inversesqrt;
#else
  {CoutLock lock; cout<<"Error: Cmatrix::inverseSqrt() called, but pMMF was compiled with no linear algebra library."<<endl;}
  exit(-1);
#endif 
}


Cvector Cmatrix::solve(const Cvector& b) const{
#ifdef _withEigen
  Eigen::MatrixXd A=convert<Eigen::MatrixXd>();
  Eigen::VectorXd be=b.convert<Eigen::VectorXd>();
  return EigenVectorXdAdaptor(A.colPivHouseholderQr().solve(be));
#elif _withLapack
	int IPIV[nrows]; 
	Cmatrix matrix_copy(*this);
	// LU factorization
	int info = LAPACKE_dgetrf(LAPACK_COL_MAJOR,nrows,ncols,matrix_copy.array,nrows,IPIV);
	if(info) {CoutLock lock; cout<<"Error: DGETRF failed!"<<endl; exit(-1);}
	Cvector x(b);
	info = LAPACKE_dgetrs(LAPACK_COL_MAJOR,'N',nrows,1,matrix_copy.array,nrows,IPIV,x.array,nrows);
	if(info) {CoutLock lock; cout<<"Error: DGETRS failed!"<<endl; exit(-1);}
	return x;
#else
  {CoutLock lock; cout<<"Error: Cmatrix::solve(...) called, but pMMF was compiled with no linear algebra library."<<endl;}
  exit(-1);
#endif 
}


pair<Cmatrix*,Cvector*> Cmatrix::symmetricEigensolver() const{
#ifdef _withEigen
  pair<Cmatrix*,Cvector*> p;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(convert<Eigen::MatrixXd>());
  p.first=new Cmatrix(EigenMatrixXdAdaptor(solver.eigenvectors()));
  p.second=new Cvector(EigenVectorXdAdaptor(solver.eigenvalues()));
  return p;
#elif _withLapack
  pair<Cmatrix*,Cvector*> p;
	p.first = new Cmatrix(*this);
	p.second = new Cvector(nrows);
	int info = LAPACKE_dsyev(LAPACK_COL_MAJOR,'V','L',nrows,p.first->array,nrows,p.second->array);
	if(info) {CoutLock lock; cout<<"Error: DSYEV failed!"<<endl; exit(-1);}
	return p;
#else
  {CoutLock lock; cout<<"Error: Cmatrix::symmetricEigensolver() called, but pMMF was compiled with no linear algebra library."<<endl;}
  exit(-1);
#endif
}


triple<Cmatrix,Cmatrix,Cvector> Cmatrix::svd() const{
#ifdef _withEigen
  triple<Cmatrix,Cmatrix,Cvector> p;
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(convert<Eigen::MatrixXd>(),Eigen::ComputeFullU|Eigen::ComputeFullV);
  p.first=Cmatrix(EigenMatrixXdAdaptor(svd.matrixU()));
  p.second=Cmatrix(EigenMatrixXdAdaptor(svd.matrixV()));
  p.third=Cvector(EigenVectorXdAdaptor(svd.singularValues()));
  return p;
#elif _withLapack

#else
  {CoutLock lock; cout<<"Error: Cmatrix::svd() called, but pMMF was compiled with no linear algebra library."<<endl;}
  exit(-1);
#endif
}


Cmatrix Cmatrix::exp(const FIELD t) const{
  return copy();
}


Cmatrix Cmatrix::eigenvectors(const int k){
#ifdef _withEigen
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(convert<Eigen::MatrixXd>());
  return Cmatrix(EigenMatrixXdAdaptor(solver.eigenvectors()),0,k);  
#else
  {CoutLock lock; cout<<"Error: Cmatrix::eigenvectors() called, but pMMF was compiled with no linear algebra library."<<endl;}
  exit(-1);
#endif
}


Cvector Cmatrix::eigenvalues(const int k) const{
#ifdef _withEigen
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(convert<Eigen::MatrixXd>());
  return Cvector(EigenVectorXdAdaptor(solver.eigenvalues()));  
#else
  {CoutLock lock; cout<<"Error: Cmatrix::eigenvectors() called, but pMMF was compiled with no linear algebra library."<<endl;}
  exit(-1);
#endif
}


// ----- In-place operations ----------------------------------------------------------------------------------

void Cmatrix::symmetrize() {
  for (int j=0;j<ncols;j++) {
    (*this)(j,j) *= 2;
    for(int i=0;i<j;i++) {
      (*this)(i,j) += (*this)(j,i);
      (*this)(j,i) = (*this)(i,j);
    }
  }
}


// ---- Iterators --------------------------------------------------------------------------------------------------


/*
CmatrixIterator Cmatrix::begin() const{ 
  return CmatrixIterator(*this,0,0,array);};

void Cmatrix::increment(CmatrixIterator& itr) const{
  itr.ptr++; itr._i++; if(itr._i>=nrows){itr._i=0; itr._j++;} }
    
CmatrixIterator Cmatrix::end() const{ 
  return CmatrixIterator(*this,0,ncols,array+nrows*ncols);}
*/


// ---- I/O -------------------------------------------------------------------------------------------------------



string Cmatrix::str(const Dense dummy) const {return Matrix::str(Dense());}
string Cmatrix::str(const Sparse dummy) const {return Matrix::str(Sparse());}
string Cmatrix::str() const {return Matrix::str(Dense());}
ostream& operator<<(ostream& stream, const Cmatrix& x){stream<<x.str(); return stream;}


string Cmatrix::classname(){return "Cmatrix";}


void Cmatrix::serialize(Rstream& rstream) const{
  rstream<<"Cmatrix{"<<Rstream::endl;
  rstream.var("nrows",nrows);
  rstream.var("ncols",ncols);
  rstream.var("array"," *OMITTED*");
  rstream<<"}"<<Rstream::endl;
}


void Cmatrix::serialize(Bofstream& ofs) const{
  ofs.tag("Cmatrix",0);
  ofs.write(nrows);
  ofs.write(ncols);
  ofs.write_array(array,nrows*ncols);
}


Cmatrix::Cmatrix(Bifstream& ifs): DenseMatrix(0,0){
  ifs.check("Cmatrix",0);
  ifs.read(nrows);
  ifs.read(ncols);
  ifs.read_array(array);
}

/*
Cmatrix::Cmatrix(DenseMatrixFile& file): Cmatrix(file.nrows,file.ncols){
  for(int i=0; i<nrows; i++)
    for(int j=0; j<ncols; j++)
      file>>array[j*nrows+i];
}


Cmatrix::Cmatrix(SparseMatrixFile& file): Cmatrix(file.nrows,file.ncols){
  for(int i=0; i<nrows*ncols; i++) array[i]=0;
  for(auto it=file.begin(); it!=file.end(); ++it) array[(*it).j*nrows+(*it).i]=(*it).value;
}
*/


Cmatrix::Cmatrix(MatrixIF& file): Cmatrix(file.nrows, file.ncols){
  file.rewind();
  if(file.sparse){
    for(int i=0; i<nrows*ncols; i++) array[i]=0;
    IndexValueTriple t;
    file>>t;
    while(t.i>=0){
      array[t.j*nrows+t.i]=t.value;
      file>>t;
    }
  }else{
    for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
	file>>array[j*nrows+i];
  }
}


void Cmatrix::saveTo(MatrixOF& file) const{
  assert(file.nrows==nrows);
  assert(file.ncols==ncols);
  if(file.sparse){
    for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
       if(array[j*nrows+i]!=0) file<<IndexValueTriple(i,j,array[j*nrows+i]);
   }else{
    for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
       file<<array[j*nrows+i];
   }
}

CSCmatrix Cmatrix::cscformat() {
  CSCmatrix result;
  result.nnz = nnz();
  result.ir = new INDEX[result.nnz];
  result.jc = new INDEX[ncols+1];
  result.val = new FIELD[result.nnz];
  result.nrows = nrows; result.ncols = ncols; 

  result.jc[ncols] = result.nnz;

  int counter=0, counter2=0;
  for(int j=0;j<ncols;j++) {
    result.jc[j]=counter;
    for (int i=0;i<nrows;i++) {
      if(array[counter2]!=0) {
        result.ir[counter] = i;
        result.val[counter] = array[counter2];
        counter++;
      }
      counter2++;
    }
  }
  return result;
} 

void Cmatrix::dump(FIELD* result) {
  std::memcpy(result,array,nrows*ncols*sizeof(FIELD));
}


// --------------------------------------------------------------------------------------------------------------



/*
Cmatrix::Cmatrix(const Cmatrix& x, const int _nrows, const int _ncols):
  Cmatrix(_nrows,_ncols){
  for(int i=0; i<nrows; i++) 
    for(int j=0; j<ncols; j++)
      array[j*nrows+i]=x(i,j);
}
*/

/*
Cmatrix(const MatrixXv& x); 
Cmatrix(const MatrixXl& x); 
Cmatrix(const MatrixXh& x); 
Cmatrix::Cmatrix(const MatrixX<Vectorv>& x): DenseMatrix(x.nrows,x.ncols){
  for(int j=0; j<ncols; j++) for(auto& p:*x.column[j]) array[j*nrows+p.first]=p.second;}
Cmatrix::Cmatrix(const MatrixX<SparseVectorl>& x): DenseMatrix(x.nrows,x.ncols){
  for(int j=0; j<ncols; j++) for(auto& p:*x.column[j]) array[j*nrows+p.first]=p.second;}
Cmatrix::Cmatrix(const MatrixX<SparseVectorh>& x): DenseMatrix(x.nrows,x.ncols){
  for(int j=0; j<ncols; j++) for(auto& p:*x.column[j]) array[j*nrows+p.first]=p.second;}
*/



