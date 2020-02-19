#ifndef SPHFD_H_
#define SPHFD_H_
#include "FDtomo/c2f.h"
#include "sphfd/vhead.h"
typedef struct{
    char timefile[160];
    char hdr[nhbyte + 1];
    float *time0;
}SPHFD_DATA;

SPHFD_DATA **sphfd(int, char**, SPEC, C2F_DATA *);
int OUTPUT_SPHFD(SPHFD_DATA *, SPEC);
#endif
