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



#include<iostream>
#include<string>
using namespace std;

class Params{
/**
A simple object that holds all the parameters necessary for the MLGkernel including
file paths for where the dataset and dataset features are stored and file path to
save the resulting gram matrix. 
**/
  public:
    // Constructer that inits the model variables.
    // set the data, feature and save paths separately.
    Params(double e, double g, int r, int l, int t, bool b): 
      eta(e), gamma(g), radius(r), levels(l), num_threads(t), grow_or_double(b) {}

  public:
    void set_paths(string data, string features){
      data_path = data; 
      features_path = features; 
    }

    void set_save_path(string save){
      save_path = save;
    }

    void show() {
      cout << "Current parameter settings:" << endl;
      cout << "    -eta            : "          << eta <<endl;
      cout << "    -gamma          : "          << gamma<<endl;
      cout << "    -radius         : "          << radius<<endl;
      cout << "    -levels         : "          << levels<<endl;
      cout << "    -grow_or_double : "          ; 
      cout << ((grow_or_double == GROW) ? "grow" : "double") <<endl;
      cout << "    -data_path      : "          << data_path<<endl;
      cout << "    -features_path  : "          << features_path<<endl;
      cout << "    -save_path      : "          << save_path<<endl;
      cout << "    -num_threads    : "          << num_threads<<endl;
    }
  public:
    static const bool GROW = 1; // subgraphs increase by given radius size at every level
    static const bool DOUBLE = 0; // subgraphs double at each level

    // MLGkernel model parameter defaults.
    double eta = 0.1;
    double gamma = 0.01;
    int radius = 1;
    int levels = 2;
    bool grow_or_double = GROW;

    // Input/Output file locations 
    string data_path;
    string features_path;
    string save_path;

    int num_threads = 1; 
};
