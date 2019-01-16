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



/**
 * To run the runMLG executable do:
 * ./runMLG -d {dataFilePath} -f {featureFilePath} -s {saveFilePath} -e {eta} -g {gamma} -r {radius} -l {levels} -t {threads}
 * Running the executable will compute the gram matrix with the MLGkernel and save it to the supplied save path.
 * Note: You must supply an absolute data/feature/save file path
 * Example usage:
 * ./runMLG -d /home/data/data.txt -f /home/data/features.txt -s /home/data/output/save.txt-e 0.1 -g 0.2 -r 3 -l 2 -t 8
**/

#include<string>
#include<sstream>
#include<stdlib.h>
#include<iostream>
#include<getopt.h>
#include "pMMFglobal.inc"
#include "MLGdataset.hpp"
#include "params.hpp"
using namespace std;

string genSaveName(string data, double eta, double gamma, int radius, int levels, bool grow_or_double) {
  stringstream ss;
  ss << fixed;
  ss.precision(2);
  ss << data << "_";

  // Round to int if needed
  if (eta == round(eta)) ss << (int)(eta+0.5) << "_";
  else ss << eta << "_";
  if (gamma == round(gamma)) ss << (int)(gamma+0.5) << "_";
  else ss << gamma << "_";

  ss << radius << "_" << levels << "_" << grow_or_double;
  ss << "_Gram";
  ss << ".txt";
  return ss.str();
}

// These are the number of discrete node labels for each of the benchmark datasets.
int get_num_features(string features){
  if(features.find("MUTAG") != string::npos) return 7;
  if(features.find("PTC") != string::npos) return 22;
  if(features.find("PROTEINS") != string::npos) return 3;
  if(features.find("NCI109") != string::npos) return 38;
  if(features.find("NCI1") != string::npos) return 37;

  cout << "Supplied dataset is not one of the sample datasets! You can manually change this code to use the correct number of discrete features of your dataset." << endl;
  exit(0);
  return 0;
}

void runMLG(Params& p) {
  threadManager.maxthreads = p.num_threads;
  MLGdataset dataset(p.data_path, p.eta, p.gamma, p.grow_or_double);
  cout << "Done loading dataset" << endl;

  if (p.features_path.empty()) {
    cout << "Computing degree features" << endl;
    for(auto g: dataset.graphs) g->computeDegreeFeatures(20); // all sample datasets have max degree < 20
  } else {
    cout << "Computing discrete features" << endl;
    int num_features = get_num_features(p.features_path);
    cout << "num features: " << num_features << endl;
    dataset.loadDiscreteFeatures(p.features_path, num_features);
  }

  dataset.computeGram(p.levels, p.radius);
  dataset.saveGram(p.save_path);
}

int main(int argc, char * argv[]){

  //threadManager.maxthreads=1;
  // Setup for parsing command line options
  const struct option longopts[] =
  {
    {"data_path", required_argument, 0, 'd'}, // path to data
    {"eta", required_argument, 0, 'e'},
    {"gamma", required_argument, 0, 'g'},
    {"radius", required_argument, 0, 'r'},
    {"levels", required_argument, 0, 'l'},
    {"data_path", required_argument, 0, 'd'}, // path to data
    {"features_path", required_argument, 0, 'f'}, // path to feature file
    {"save_path", required_argument, 0, 's'}, // location to save kernel file
    {"threads", required_argument, 0, 't'}, // number of threads to use
    {"grow_or_double", required_argument, 0, 'm'}, // bool flag to grow(by the radius) or double
    {0,0,0,0}
  };

  // Some default settings
  int index;
  int iarg=0;
  string data_path, save_path, features_path;
  double eta = 0.1;
  double gamma = 0.01;
  int radius = 2;
  int levels = 2;
  int threads = 1; // note the max threads you can use depends on your machine
  bool grow_or_double = Params::GROW;
  while(iarg != -1)
  {
    iarg = getopt_long(argc, argv, "d:e:g:r:l:f:t:m:s:", longopts, &index);
    switch (iarg) {
      case 'd':
        data_path = optarg;
        break;
      case 'e':
        eta = atof(optarg);
        break;
      case 'f':
        features_path = optarg;
        break;
      case 'g':
        gamma = atof(optarg);
        break;
      case 'r':
        radius = atoi(optarg);
        break;
      case 'l':
        levels = atoi(optarg);
        break;
      case 's':
        save_path = optarg;
        break;
      case 'm':
        grow_or_double = atoi(optarg); // 1 to grow by leaf radius, 0 to double
        break;
      case 't':
        threads = atoi(optarg);
        break;
    }
  }

  Params params(eta, gamma, radius, levels, threads, grow_or_double);
  if (data_path.empty()) cout << "You must supply a data file to run MLG." << endl;
  if (save_path.empty()) {
    save_path = genSaveName(data_path, eta, gamma, radius, levels, grow_or_double);
    cout << "Save path not provided, so we will save the gram matrix to: " <<  save_path << endl;
  }
  params.set_paths(data_path, features_path);
  params.set_save_path(save_path);
  params.show();

  runMLG(params);
  cout << "Done with mlg kernel" << endl;
}
