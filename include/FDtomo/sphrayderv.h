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
    int total_elements;
}sparse_matrix;

typedef struct{
    sparse_matrix *mat;
    float *b;
}SPHRAYDERV_DATA;

typedef struct{
    int iray, iraystat, idatout, nomat, dmean, kmin, kmax, ido1d;
    float vpvsscale, resflag;
    //msut files
    char locdfil[MAXSTRLEN + 1], telrerr[MAXSTRLEN + 1], 
    	dtdsfil[MAXSTRLEN + 1], resfile[MAXSTRLEN + 1], hitfile[MAXSTRLEN + 1], 
		dtdhfil[MAXSTRLEN + 1], bookfil[MAXSTRLEN + 1], sclefil[MAXSTRLEN + 1];
    //optinal files
    char telefil[MAXSTRLEN],  pbasfil[MAXSTRLEN], sbasfil[MAXSTRLEN], shotfil[MAXSTRLEN],
		elipfil[MAXSTRLEN], raystat[MAXSTRLEN],  dotfile[MAXSTRLEN], 
		headfil[MAXSTRLEN], entfile[MAXSTRLEN], stcfile[MAXSTRLEN], specfile[MAXSTRLEN];
}SphraydervEnv;

SPHRAYDERV_DATA *sphrayderv(velocityModel3D, travelTimeTable *, Event *, int, Station *, int, SphraydervEnv, CommonEnv);
int OUTPUT_SPHFRAYDERV(SPHRAYDERV_DATA *, SPEC);
int LOG_SPHFRAYDERV(SPEC);
SphraydervEnv setSphraydervEnv(char *);
void setSphrayderVariables(SphraydervEnv *, char *);
void setSphrayderFiles(SphraydervEnv *, char *);
#endif
