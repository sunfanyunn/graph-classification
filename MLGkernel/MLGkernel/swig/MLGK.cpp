#include "ThreadManager.hpp"
#include "pMMFbase.hpp"
#include "Vectorv.hpp"
#include "Vectorl.hpp"
#include "Vectorh.hpp"

std::default_random_engine randomNumberGenerator;
bool multithreading=true;
ThreadManager threadManager(4);
char strbuffer[255];
mutex CoutLock::mx;
