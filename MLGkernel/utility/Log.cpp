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

#include "Log.hpp"
#include <iostream>
#include <mutex>

//extern mutex cout_mutex;


Log& Log::operator<<(const string& s){
  CoutLock lock;
  if(stream==nullptr) cout<<s<<endl;
  else (*stream)<<s<<"\n";
  skippedlines=0;
  return *this; 
}

Log& Log::operator<<(const char* s){
  CoutLock lock;
  double now=chrono::duration<double>(chrono::system_clock::now()-t).count();  
  if(stream==nullptr){
    cout.precision(5);
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout<<now<<"s: "<<s<<endl;
  }else{
    (*stream)<<s<<"\n";
  }
  skippedlines=0;
  return *this; 
}


Log& Log::log(const int v, const string& s){
  if(verbosity<v) return *this;
  CoutLock lock;
  double now=chrono::duration<double>(chrono::system_clock::now()-t).count(); 
  if(stream==nullptr){
    cout.precision(5);
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout<<now<<"s: "<<s<<endl;
  }else{
    (*stream)<<now<<"s: "<<s<<"\n";
  }
  skippedlines=0;
  return *this;
}


Log& Log::log(const int v, const char* s, const int i){
  if(verbosity<v) return *this;
  CoutLock lock;
  double now=chrono::duration<double>(chrono::system_clock::now()-t).count();  
  if(stream==nullptr){
    cout.precision(5);
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout<<now<<"s: "<<s<<i<<endl;
  }else{
    (*stream)<<now<<"s: "<<s<<"\n";
  }
  skippedlines=0;
  return *this;
}


Log& Log::log(const int v, const char* s1, const int i, const char* s2){
  if(verbosity<v) return *this;
  CoutLock lock;
  double now=chrono::duration<double>(chrono::system_clock::now()-t).count();  
  if(stream==nullptr){
    cout.precision(5);
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout<<now<<"s: "<<s1<<i<<s2<<endl;
  }else{
    (*stream)<<now<<"s: "<<s1<<i<<s2<<"\n";
  }
  skippedlines=0;
  return *this;
}


Log& Log::log(const int v, const char* s1, const int i1, const char* s2, const int i2, const char* s3){
  if(verbosity<v) return *this;
  CoutLock lock;
  double now=chrono::duration<double>(chrono::system_clock::now()-t).count();  
  if(stream==nullptr){
    cout.precision(5);
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout<<now<<"s: "<<s1<<i1<<s2<<i2<<s3<<endl;
  }else{
    (*stream)<<now<<"s: "<<s1<<i1<<s2<<i2<<s3<<"\n";
  }
  skippedlines=0;
  return *this;
}


Log& Log::log(const int v, const char* s1, const int i1, const char* s2, const int i2, const char* s3, const int i3, const char* s4){
  if(verbosity<v) return *this;
  CoutLock lock;
  double now=chrono::duration<double>(chrono::system_clock::now()-t).count();
  if(stream==nullptr){
    cout.precision(5);
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout<<now<<"s: "<<s1<<i1<<s2<<i2<<s3<<i3<<s4<<endl;
  }else{
    (*stream)<<now<<"s: "<<s1<<i1<<s2<<i2<<s3<<i3<<s4<<"\n";
  }
  skippedlines=0;
  return *this;
}


Log& Log::log(const int v, const char* s1, const int i1, const char* s2, const int i2, const char* s3, const int i3, const char* s4, const int i4, const char* s5){
  if(verbosity<v) return *this;
  CoutLock lock;
  double now=chrono::duration<double>(chrono::system_clock::now()-t).count();  
  if(stream==nullptr){
    cout.precision(5);
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout<<now<<"s: "<<s1<<i1<<s2<<i2<<s3<<i3<<s4<<i4<<s5<<endl;
  }else{
    (*stream)<<now<<"s: "<<s1<<i1<<s2<<i2<<s3<<i3<<s4<<i4<<s5<<"\n";
  }
  skippedlines=0;
  return *this;
}


Log& Log::log(const int v, const char* s, const double f){
  if(verbosity<v) return *this;
  CoutLock lock;
  double now=chrono::duration<double>(chrono::system_clock::now()-t).count();  
  if(stream==nullptr){
    cout.precision(5);
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout<<now<<"s: "<<s<<f<<endl;
  }else{
    (*stream)<<now<<"s: "<<s<<f<<"\n";
  }
  skippedlines=0;
  return *this;
}

Log& Log::log(const int v, const int i, const double f){
  if(verbosity<v) return *this;
  CoutLock lock;
  double now=chrono::duration<double>(chrono::system_clock::now()-t).count();  
  if(stream==nullptr){
    cout.precision(5);
    cout.setf(ios_base::fixed, ios_base::floatfield);
    cout<<now<<" "<<i<<" "<<f<<" "<<endl;
  }else{
    (*stream)<<now<<" "<<i<<" "<<f<<" "<<"\n";
  }
  skippedlines=0;
  return *this;
}

void Log::startClock(const int i){
  //while(i>time.size()-1) time.push_back(chrono::system_clock::now());
  t=chrono::system_clock::now();
  //*this<<"Start clock";
  //time[i]=new chrono::time_point<chrono::system_clock>(chrono::system_clock::now());
}

double Log::clock(const int i){
  return chrono::duration<double>(chrono::system_clock::now()-t).count();
}
