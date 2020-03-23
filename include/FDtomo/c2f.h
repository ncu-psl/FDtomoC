#ifndef C2F_H_
#define C2F_H_
#include "common/read_spec.h"
#include "FDtomo/make1d.h"
typedef struct{
    char filename[MAXSTRLEN + 1];
    struct vhead head;    
    float vsave[nxyzm2];
}VELFILE;

typedef struct{
    VELFILE *vpfile;
    VELFILE *vsfile;
}C2F_DATA;

C2F_DATA *c2f(SPEC, MAKE1D_DATA *);
int OUTPUT_C2F(C2F_DATA *, SPEC);
#endif
