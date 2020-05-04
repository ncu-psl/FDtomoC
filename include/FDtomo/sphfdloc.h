#ifndef SPHFDLOC_H_
#define SPHFDLOC_H_
#include "common/earthquake.h"
#include "common/travel_time.h"
#include "FDtomo/sphfd.h"
typedef struct{
    char event_hdr[100];
    char event[20000];
}SPHFDLOC_DATA;

LocData *sphfdloc(Coordinate3D, travelTimeTable *, EventNode *, StationNode *, LocEnv);
int OUTPUT_SPHFDLOC(SPHFDLOC_DATA **, SPEC);
int LOG_SPHFDLOC(SPEC);
#endif
