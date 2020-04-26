#ifndef INTERPOLATION
#define INTERPOLATION
#include "common/read_spec.h"
#include "common/grid.h"
float linear_interpolation(float, float, float, float, float);
float *linear_interpolation_array(float *, float *, float*, int, int, char *);
float trilinear_interpolation(Point3D, Cell cells[2][2][2]);
#endif