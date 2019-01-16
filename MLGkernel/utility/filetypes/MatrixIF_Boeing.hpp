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
#ifndef _MatrixIF_Boeing
#define _MatrixIF_Boeing

#include "MatrixIF.hpp"
//#include "hb_io.hpp"


class MatrixIF_Boeing: public MatrixIF{
public:

  MatrixIF_Boeing(const string filename){
    sparse=1;
    ifs.open(filename);
    if(ifs.fail()){cout<<"Failed to open "<<filename<<"."<<endl; return;}

    int total_num_lines,num_pointer_lines, num_rowidx_lines, num_value_lines, num_rhs_lines=0, dummy;
    char matrix_type[4];
    bool rhs_header=0;
    ifs.ignore(numeric_limits<streamsize>::max(),'\n'); //skip to second line
    ifs>>total_num_lines>>num_pointer_lines>>num_rowidx_lines>>num_value_lines;
    if (ifs.peek()!='\n') ifs>>num_rhs_lines;
    ifs>>matrix_type>>nrows>>ncols>>nnz;
    
    if(ifs.fail()) {cout<<"Corrupt header in file!"<<endl; return;}
    
    if(strchr("RrPp",matrix_type[0])==NULL || strchr("Ss",matrix_type[1])==NULL || strchr("Aa",matrix_type[2])==NULL) {
      cout<<"Matrix should be real-valued, symmetric and sparse!"<<endl; 
      return;
    }
    ifs.ignore(numeric_limits<streamsize>::max(),'\n'); //skip to next line
    ifs.ignore(numeric_limits<streamsize>::max(),'\n'); //skip info about formats
    if(num_rhs_lines) {rhs_header=1; ifs.ignore(numeric_limits<streamsize>::max(),'\n');} //ignore header info about right-hand side vectors

    rowfs.open(filename);
    colfs.open(filename);
    valfs.open(filename);
    for(int i=0;i<4+rhs_header && colfs.good(); i++) 
      colfs.ignore(numeric_limits<streamsize>::max(),'\n'); 
    for(int i=0;i<4+rhs_header+num_pointer_lines && rowfs.good(); i++) 
      rowfs.ignore(numeric_limits<streamsize>::max(),'\n'); 
    for(int i=0;i<4+rhs_header+num_pointer_lines+num_rowidx_lines && valfs.good(); i++) 
      valfs.ignore(numeric_limits<streamsize>::max(),'\n'); 

    if (!colfs.good() || !rowfs.good() || !valfs.good()) { cout<<"Corrupt file (header part is ok)!"<<endl; return; }

    colg = colfs.tellg();
    rowg = rowfs.tellg();
    valg = valfs.tellg();

    colfs>>dummy>>colptr_value;
    nnz_in_col = colptr_value-dummy;
    rowcounter=0;
    colcounter=1;
    nnzcounter=0;

    secondSweep=false;

    if (strchr("Pp",matrix_type[0])!=NULL) pointers_only=true; else pointers_only=false;
  }


  ~MatrixIF_Boeing(){ rowfs.close(); colfs.close(); valfs.close(); }

public:

  MatrixIF& operator>>(IndexValueTriple& dest){
    if(!rowfs.good() || !colfs.good() || !valfs.good()) {
      cout<<"end-of-file reached."<<endl;
      return *this;
    }
    if(nnzcounter>=nnz) {
      if (!secondSweep) {
        secondSweep=true;
        rowcounter=0;
        colcounter=1;
        nnzcounter=0;
        colfs.seekg(colg);
        rowfs.seekg(rowg);
        valfs.seekg(valg);
        int dummy;
        colfs>>dummy>>colptr_value;
        nnz_in_col = colptr_value-dummy;
      } else dest.i=-1;
      return *this;
    }
    if(rowcounter>=nnz_in_col) {
      int dummy = colptr_value;
      colfs>>colptr_value;
      nnz_in_col = colptr_value-dummy;
      colcounter++;
      rowcounter=0;
    }
    if (nnz_in_col > 0) {
      rowfs>>dest.i;
      dest.j=colcounter;
      if (pointers_only) dest.value=1;
      else valfs>>dest.value;
      dest.i--; dest.j--;
      if (secondSweep) { int dummy; dummy=dest.i; dest.i=dest.j; dest.j=dummy;}
      rowcounter++;
      nnzcounter++;
    }
    return *this;
  }


public:
  ifstream rowfs,colfs,valfs;
  streampos rowg,colg,valg;
  int nnz_in_col, rowcounter, colcounter;
  int nnz, nnzcounter;
  int colptr_value;
  bool pointers_only;
  bool secondSweep;
};


#endif
