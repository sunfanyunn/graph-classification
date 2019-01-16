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


#include "pMMFbase.hpp" 
#include "Vectorv.hpp"
#include "Vectorl.hpp"
#include "Vectorh.hpp"
#include "ThreadManager.hpp"

bool multithreading=true;
ThreadManager threadManager(4);
mutex cout_mutex;
mutex CoutLock::mx;

std::default_random_engine randomNumberGenerator;

FIELD Vectorv::dummyZero=0;
FIELD Vectorl::dummyZero=0;
FIELD Vectorh::dummyZero=0;

char strbuffer[255];

Log mlog;

