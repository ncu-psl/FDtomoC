#ifndef VELOCITY_MODEL
#define VELOCITY_MODEL
#include "read_spec.h"
#include "common/vec/vec.h"
#include "common/grid.h"
#include "common/parameter.h"
#include "common/interpolation.h"
#include "FDtomo/make1d.h"

#define MAX1D 1000
typedef struct velocityModel1D_{
    Coordinate1D coordinate;
    float *velocity;
}velocityModel1D;

typedef struct velocityModel3D_{
    Coordinate3D coordinate;
    float *velocity;
}velocityModel3D;

void readVelocityModel1D(char *, velocityModel1D *, velocityModel1D *, char *);
void transform1D(Coordinate1D, velocityModel1D *,char *);
velocityModel3D create3DModel(Coordinate3D, velocityModel1D);
velocityModel3D transform3D(Coordinate3D, velocityModel3D *);
float getPointVel(Point3D, velocityModel3D *);
Point3D getPoint3DModel(Point3D, velocityModel3D *);
float trilinear_interpolation_base(Point3D , Point3D, Coordinate3D, velocityModel3D *);
velocityModel3D change2ColumnMajor(velocityModel3D);
void output3DModel(velocityModel3D, char *);

#endif