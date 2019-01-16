/*
 -----------------------------------------------------------------------------
 
 MLGkernel is an open source implementation of the Multiscale Laplacian Graph
 Kernel for computing the gram matrix of a collection of graphs.
 
 Copyright (C) 2016 Imre Risi Kondor, Horace Pan
 
 
 This program is free software; you can redistribute it and/or
 modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2
 of the License, or (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, see <http://www.gnu.org/licenses/>.
 
 ----------------------------------------------------------------------------- */



#ifndef _MLGdataset
#define _MLGdataset

#include "MLGgraph.hpp"
#include <string>

class MLGdataset{
public:

  MLGdataset(){}
  MLGdataset(const std::string filename, double eta, double gamma, bool grow): gamma(gamma), grow(grow), eta(eta){
    loadGraphs(filename);
  }
  ~MLGdataset() {for(auto p:graphs) delete p;}

public:

  void condense(const int nlevels, const int leaf_radius=2);
  void computeGram(const int levels, const int radius);

public:

  void loadGraphs(std::string filename);
  void loadDiscreteFeatures(std::string filename, int numFeatures);
  void loadFeatures(std::string filename);
  void saveGram(std::string filename);
  void fillGram(double *npmatrix, int rows, int cols);

public:

  vector<MLGgraph*> graphs;
  double gamma; // regularizer constant
  double eta; // regularizer constant
  int levels;
  int radius;
  bool grow; // 1 to grow by the leaf radius, 0 to double
  Cmatrix gram;
};

#endif
