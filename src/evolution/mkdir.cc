#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>

// function to create folder for LHAPDFgrid.f
extern "C" void mkdir_(char* name, int) { mkdir(name, 0777); }
