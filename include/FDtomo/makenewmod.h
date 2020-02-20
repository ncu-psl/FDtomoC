#ifndef MAKENEWMOD_H_
#define MAKENEWMOD_H_
#include "FDtomo/runlsqr.h"
#define nhbyte 58 * 4

typedef struct{
    char hdr[nhbyte + 1];
    int igridx[nxcm1], igridy[nycm1], igridz[nzcm1];
    float *vn;

}MAKENEWMOD_DATA;
MAKENEWMOD_DATA *makenewmod(SPEC, RUNLSQR_DATA *);
int OUTPUT_MAKENEWMOD(MAKENEWMOD_DATA *, SPEC);

#endif
