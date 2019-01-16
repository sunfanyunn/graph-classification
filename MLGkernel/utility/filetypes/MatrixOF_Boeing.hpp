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
#ifndef _MatrixOF_Boeing
#define _MatrixOF_Boeing

#include "MatrixOF.hpp"
#include <iomanip>

class MatrixOF_Boeing: public MatrixOF{
public:

  MatrixOF_Boeing(const string filename, const int _nrows, const int _ncols){
    nrows=_nrows; ncols=_ncols; sparse=1; 
    
    // Write row indices, column pointers, values in separate files and concatenate
    rowfs.open("row.temp");
    colfs.open("col.temp");
    valfs.open("val.temp");

    nnz=0;
    nnz_in_col=0;
    rowchars=colchars=valchars=0;
    colptr_val=1;
    
    colfs<<colptr_val;
    
    current_col=0;
    num_rowidx_lines=num_colptr_lines=num_value_lines=1;
    filename_local=filename;
  }
  
  ~MatrixOF_Boeing() { 
    rowfs<<endl;
    colfs<<" "<<colptr_val+nnz_in_col<<endl;
    valfs<<endl;
    rowfs.close();
    colfs.close();
    valfs.close();

    // write header
    ofs.open(filename_local);
    ofs<<"Matrix"<<setw(72)<<"RSA_32"<<endl;
    ofs<<setw(14)<<num_colptr_lines+num_rowidx_lines+num_value_lines<<setw(14)<<num_colptr_lines<<setw(14)<<num_rowidx_lines<<setw(14)<<num_value_lines<<setw(14)<<0<<endl;
    ofs<<"RSA"<<setw(25)<<nrows<<setw(14)<<ncols<<setw(14)<<nnz<<setw(14)<<0<<endl;
    ofs<<"(16I5)          (16I5)          (10F7.1)            (10F7.1)"<<endl;

    // write data
    ifstream tempfs;
    tempfs.open("col.temp");
    ofs<<tempfs.rdbuf();
    tempfs.close(); tempfs.open("row.temp");
    ofs<<tempfs.rdbuf();
    tempfs.close(); tempfs.open("val.temp");
    ofs<<tempfs.rdbuf();
    tempfs.close();

    // delete temporary files
    remove("row.temp");
    remove("col.temp");
    remove("val.temp");
  }

public:

  MatrixOF& operator<<(const IndexValueTriple& t){
    if (t.i>t.j) return *this; // write only the upper triangular part
    if(t.j==current_col+1) { // no more nonzeros in previous column
      if (colchars+1+std::to_string(colptr_val+nnz_in_col).length() > 80) {
	colfs<<endl;
	colchars=0;
	num_colptr_lines++;
      }
      current_col++;
      colfs<<" "<<colptr_val+nnz_in_col;
      colptr_val += nnz_in_col;
      colchars += 1+std::to_string(colptr_val+nnz_in_col).length();
      nnz_in_col=0;
    }
    if (t.j==current_col) {
      nnz_in_col++;
      nnz++;
      if(rowchars+1+std::to_string(t.i+1).length() > 80) {
	rowfs<<endl;
	rowchars=0;
	num_rowidx_lines++;
      }
      rowfs<<" "<<t.i+1;
      rowchars += 1+std::to_string(t.i+1).length();
      if(valchars+1+std::to_string(t.value).length() > 80) {
	valfs<<endl;
	valchars=0;
	num_value_lines++;
      }
      valfs<<" "<<t.value;
      valchars += 1+std::to_string(t.value).length();
    } else { cout <<"Matrix should be written column-wise!!"<<current_col<<" "<<t.j<<endl; return *this; }
    
    return *this;
  }

public:
  ofstream headerfs,rowfs,colfs,valfs;
  int nnz_in_col,nnz,colptr_val,current_col;
  int rowchars,colchars,valchars;
  int num_rowidx_lines,num_colptr_lines,num_value_lines;
  string filename_local;

};



#endif
