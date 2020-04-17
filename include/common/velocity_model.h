#ifndef VELOCITY_MODEL
#define VELOCITY_MODEL
#include "read_spec.h"
#include "common/grid.h"
#include "FDtomo/make1d.h"

#define MAX1D 1000

typedef struct{
    float vp[MAX1D], vs[MAX1D], z[MAX1D];
    char terp[MAX1D + 1];
	int nl;
}velocity1D;

typedef struct{
    Mesh mesh;
    float *vp, *vs;
}velocity3D;

velocity1D read_velocity1D(SPEC spec);
velocity3D create3DModel(Mesh, velocity1D);
velocity3D generate3DModel(float *, float *, Mesh);
velocity3D transform(velocity3D);
float getPointVel(Point3D, velocity3D, char);
#endif