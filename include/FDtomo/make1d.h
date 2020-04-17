#ifndef MAKE1D_H_
#define MAKE1D_H_
#include "common/read_spec.h"
#include "common/vhead.h"
#include "common/parameter.h"
#define nhbyte 58 * 4

typedef struct{
    struct vhead head;
    int igridx[nxcm1], igridy[nycm1], igridz[nzcm1];
    float vsave[nxyzcm2];

}MAKE1D_DATA;


float flatvel(float, float);
float uflatz(float);
float flatz(float);
char * dtoa(char *, double, int);
int OUTPUT_MAKE1D(MAKE1D_DATA *, SPEC);
int LOG_MAKE1D(SPEC);

#endif
