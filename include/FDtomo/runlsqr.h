#ifndef RUNLSQR_H_
#define RUNLSQR_H_
#include "FDtomo/sphrayderv.h"
typedef struct{
    int n;
    float *x;
    int *jndx;
    float *se;
}RUNLSQR_DATA;

typedef struct{
    int intlim;
    float damper;
    char nmodfil[MAXSTRLEN + 1], fresfil[MAXSTRLEN + 1];
}RunlsqrEnv;


RUNLSQR_DATA *runlsqr(SPHRAYDERV_DATA *, RunlsqrEnv, CommonEnv);
int OUTPUT_RUNLSQR(RUNLSQR_DATA *, SPEC);
int LOG_RUNLSQR(SPEC);
RunlsqrEnv setRunlsqrEnv(char *);
void setRunlsqrVariables(RunlsqrEnv *,char *);
void setRunlsqrFiles(RunlsqrEnv *, char *);


#endif
