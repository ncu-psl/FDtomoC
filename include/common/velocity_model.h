#ifndef VELOCITY_MODEL
#define VELOCITY_MODEL
#include "read_spec.h"
#include "FDtomo/make1d.h"

#define MAX1D 1000

typedef struct{
    float vp[MAX1D], vs[MAX1D], z[MAX1D];
    char terp[MAX1D + 1];
	int nl;
}velocity1D;

typedef struct{
    GRID grid;
    float *vp, *vs;
}velocity3D;

velocity1D read_velocity1D(SPEC spec);
velocity3D create3DModel(SPEC, velocity1D);
velocity3D generate3DModel(float *, float *, GRID);

#endif