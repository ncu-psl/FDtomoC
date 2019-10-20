#ifndef SHARED_VARIABLES
#define SHARED_VARIABLES
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
#endif