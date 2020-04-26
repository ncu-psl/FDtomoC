#ifndef SPHFDLOC_H_
#define SPHFDLOC_H_
#include "common/earthquake.h"
#include "FDtomo/sphfd.h"
typedef struct{
    char event_hdr[100];
    char event[20000];
}SPHFDLOC_DATA;

//SPHFDLOC_DATA **sphfdloc(Mesh, Time, Event);
int OUTPUT_SPHFDLOC(SPHFDLOC_DATA **, SPEC);
int LOG_SPHFDLOC(SPEC);
#endif
