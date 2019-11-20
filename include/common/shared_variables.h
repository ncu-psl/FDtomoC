#ifndef SHARED_VARIABLES
#define SHARED_VARIABLES
#define MAXSTRLEN 132
extern char VERSION[10];
extern double rearth, degrad, hpi;
extern char timedir[60 + 1];
extern char eqkdir[60 + 1];
extern int ittnum;
extern double clat, clon, cz;
extern float az;
extern int iflat, isph, vs1d;
extern float azmod;

//sphfdloc
extern int iread, ivs, nthres, kmin;
extern double vpvs, resthres, resthrep, stdmax;
extern int ndiv, ndiv2, total_earthquakes;
extern int num_parfiles;

//sphrayderv
extern int ivpvs, idmean, iray, iraystat, idatout, nomat;
extern int istacor, idoshot, idotel, kmax;
extern float vpvsscale, resflag;

//runlsqr
extern int intlims;
extern float damper;

//makenewmod
extern int limitu, ivpvs, mavx, mavy, mavz;
extern int nsmooth, ipscflg, ido1d;
extern float dvperc, pertscl;

//files
extern char oldvfil[MAXSTRLEN + 1], onedfil[MAXSTRLEN + 1]; //make1d
extern char tgrdfil[MAXSTRLEN + 1], finevel[MAXSTRLEN + 1]; //c2f
extern char leqsfil[MAXSTRLEN + 1], fsumfil[MAXSTRLEN + 1], outlfil[MAXSTRLEN + 1],
             fhedfil[MAXSTRLEN + 1], fdatfil[MAXSTRLEN + 1]; //sphfdloc
extern char stafile[MAXSTRLEN + 1], locdfil[MAXSTRLEN + 1], telrerr[MAXSTRLEN + 1], 
            dtdsfil[MAXSTRLEN + 1], resfile[MAXSTRLEN + 1], hitfile[MAXSTRLEN + 1], 
            dtdhfil[MAXSTRLEN + 1], bookfil[MAXSTRLEN + 1], sclefil[MAXSTRLEN + 1];  //sphrayderv
extern char nmodfil[MAXSTRLEN + 1], fresfil[MAXSTRLEN + 1]; // runlsqr
extern char fmodfil[MAXSTRLEN + 1];

#endif