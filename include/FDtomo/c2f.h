#ifndef C2F_H_
#define C2F_H_
#include "common/read_spec.h"
#include "FDtomo/make1d.h"
typedef struct{
    char hdr[nhbyte + 1];
    float vsave[nxyzm2];
}VELFILE;

typedef struct{
    VELFILE *vpfile;
    VELFILE *vsfile;
}C2F_DATA;

C2F_DATA *c2f(SPEC, make1d_data *);

#endif
