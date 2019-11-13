#ifndef MAKE1D_H_
#define MAKE1D_H_
#include "common/read_spec.h"
#include "parameter.h"
#define nhbyte 58 * 4

typedef struct{
    char hdr[nhbyte + 1];
    int igridx[nxcm1], igridy[nycm1], igridz[nzcm1];
    float vsave[nxyzcm2];

}make1d_data;
make1d_data *make1d(SPEC);

#endif
