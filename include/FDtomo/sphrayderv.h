#ifndef SPHRAYDERV_H_
#define SPHRAYDERV_H_
#include "FDtomo/sphfdloc.h"
typedef struct{
    float elements[SIZEOFA];
    int column_elements[SIZEOFA];
    int elements_row[NMAX];
    int *jndx;
    int number_columns;
    int number_rows;
}sparse_matrix;

typedef struct{
    sparse_matrix *mat;
    float *b;
}SPHRAYDERV_DATA;

SPHRAYDERV_DATA *sphrayderv(SPEC, SPHFDLOC_DATA **);

#endif
