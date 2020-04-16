#ifndef READ_SPEC_H
#define READ_SPEC_H
#include "parseprogs.h"
#include "common/shared_variables.h"
#include "common/gridspec.h"
#define MAX1D 1000
#define MAXSTRLEN 132
typedef struct 
{
char spec_file[MAXSTRLEN];

GRID grid;

double clat, clon, cz;
float az, azmod;
int iflat, isph, vs1d;

//sphfdloc
int iread;
int ivs;
double vpvs;
int nthres;
double resthres;
double resthrep;
double stdmax;
int kmin;
int ndiv;
int ndiv2;
int ittnum;
int total_earthquakes;
char timedir[60 + 1];
char eqkdir[60 + 1];

//sphrayderv
float vpvsscale;
int idmean, iray, iraystat, idatout, nomat;
int ivpvs, istacor, idoshot, idotel, kmax;
float resflag;

//runlsqr
float damper;
int intlims;

//makenewmod
int limitu, mavx, mavy, mavz;
int nsmooth, ipscflg, ido1d;
float dvperc, pertscl;

//must files
char oldvfil[MAXSTRLEN + 1], onedfil[MAXSTRLEN + 1]; // make1d
char tgrdfil[MAXSTRLEN + 1], finevel[MAXSTRLEN + 1]; // c2f
char leqsfil[MAXSTRLEN + 1], fsumfil[MAXSTRLEN + 1], outlfil[MAXSTRLEN + 1],
             fhedfil[MAXSTRLEN + 1], fdatfil[MAXSTRLEN + 1]; //sphfdloc
char stafile[MAXSTRLEN + 1], locdfil[MAXSTRLEN + 1], telrerr[MAXSTRLEN + 1], 
    dtdsfil[MAXSTRLEN + 1], resfile[MAXSTRLEN + 1], hitfile[MAXSTRLEN + 1], 
	dtdhfil[MAXSTRLEN + 1], bookfil[MAXSTRLEN + 1], sclefil[MAXSTRLEN + 1];  //sphrayderv
char nmodfil[MAXSTRLEN + 1], fresfil[MAXSTRLEN + 1]; //runlsqr			
char fmodfil[MAXSTRLEN + 1]; //makenewmod
char parlist[MAXSTRLEN + 1];

//optional files parameter
//sphrayderv
char telefil[MAXSTRLEN],  pbasfil[MAXSTRLEN], sbasfil[MAXSTRLEN], shotfil[MAXSTRLEN],
		elipfil[MAXSTRLEN], raystat[MAXSTRLEN],  dotfile[MAXSTRLEN], 
		headfil[MAXSTRLEN], entfile[MAXSTRLEN], stcfile[MAXSTRLEN], specfile[MAXSTRLEN];
// makenewmod
char nstafil[MAXSTRLEN + 1];

}SPEC;

typedef struct{
    float vp[MAX1D][2], z[MAX1D];
    char terp[MAX1D + 1];
	int nl;
}velocity1D;

void read_variables(char *spec_file, SPEC *spec);
void read_files(char *spec_file, SPEC *file_identifier);
void read_grid(char *spec_file, SPEC *spec);
void read_error(char *name, char *type, FILE *fp_spc);
velocity1D read_velocity1D(SPEC spec);
SPEC create_spec(char *specFile);
#endif
