#ifndef RUNLSQR_H_
#define RUNLSQR_H_
#include "FDtomo/sphrayderv.h"
typedef struct{
    int n;
    float *x;
    int *jndx;
    float *se;
}RUNLSQR_DATA;

RUNLSQR_DATA *runlsqr(SPEC, SPHRAYDERV_DATA *);

#endif
