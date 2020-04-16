#ifndef MAKE1D_H_
#define MAKE1D_H_
#include "common/read_spec.h"
#include "common/vhead.h"
#include "parameter.h"
#define nhbyte 58 * 4

typedef struct{
    struct vhead head;
    int igridx[nxcm1], igridy[nycm1], igridz[nzcm1];
    float vsave[nxyzcm2];

}MAKE1D_DATA;

typedef struct{
    GRID grid;
    float *vp, *vs;
}velocity3D;


velocity3D create3DModel(SPEC, velocity1D);

velocity3D generate3DModel(float *, float *, GRID);
int OUTPUT_MAKE1D(MAKE1D_DATA *, SPEC);
int LOG_MAKE1D(SPEC);

#endif
