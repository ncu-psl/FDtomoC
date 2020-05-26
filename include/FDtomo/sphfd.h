#ifndef SPHFD_H_
#define SPHFD_H_
#include "FDtomo/c2f.h"
#include "common/grid.h"
#include "common/station.h"
#include "common/velocity_model.h"
#include "common/travel_time.h"
#include "common/vhead.h"
typedef struct{
    char timefile[160];
    char hdr[nhbyte + 1];
    float *time0;
}SPHFD_DATA;

travelTimeTable *sphfd(velocityModel3D, Station *, int);
travelTimeTable sphfd_exec(velocityModel3D, Station);
int OUTPUT_SPHFD(SPHFD_DATA *, SPEC);
#endif
