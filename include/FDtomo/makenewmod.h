#ifndef MAKENEWMOD_H_
#define MAKENEWMOD_H_
#include "FDtomo/runlsqr.h"
#include "common/velocity_model.h"
#define nhbyte 58 * 4

typedef struct{
    struct vhead head;
    int igridx[nxcm1], igridy[nycm1], igridz[nzcm1];
    float *vn;

}MAKENEWMOD_DATA;

typedef struct{
    int mavx, mavy, mavz, nsmooth, limitu, ipscflg, ido1d;
    float dvperc, pertscl;

}MakenewmodEnv;

MAKENEWMOD_DATA *makenewmod(Coordinate3D, velocityModel3D, velocityModel3D, RUNLSQR_DATA *, int, MakenewmodEnv, CommonEnv);
int OUTPUT_MAKENEWMOD(MAKENEWMOD_DATA *, SPEC);
int LOG_MAKENEWMOD(SPEC);
MakenewmodEnv setMakeNewmodEnv(char *);
void setMakeNewmodVariables(MakenewmodEnv *, char *);
#endif
