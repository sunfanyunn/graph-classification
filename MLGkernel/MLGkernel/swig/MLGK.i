%module MLGK
%include "std_string.i"
%{
  #define SWIG_FILE_WITH_INIT
  #include "../params.hpp"
  #include "../MLGdataset.hpp"
  using namespace std;
%}

%include "numpy.i"
%init
%{
import_array();
%}

%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *npmatrix, int rows, int cols)};
%include "../params.hpp"
%include "../MLGdataset.hpp"
