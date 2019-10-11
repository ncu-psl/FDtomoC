// c	Arrival time inversion code for use in conjunction with FD travel
// c	time tables.This is a version of telrayderv.f that uses spherical
// c	coordinates.
// c
// c	This program accumulates derivatives in a "metagrid"
// c	in which a coarse grid is overlain on the original fine grid.
// c
// c	In this version the fine and coarse grid points overlap exactly, so we specify the
// c	the coarse grid in terms of multiples of the fine grid spacing "h".In this version
// c	the grid spacing in any direction(hx, hy, hz) is a function ONLY of that particular
// c	direction(hx(x), hy(y), hz(z)).There is no reason why these can't be, but it's easier
// c	bookkeeping if they aren't.
// c
// c	We specify combinations of grid points by the number of CELLS they subsume.For example,
// c	suppose we have a fine grid with 24 points in the x direction, which defines 23 cells.
// c	A combination of cells could be(6, 4, 2, 1, 4, 6) - the sum must be 23 - which results in
// c	7 grid points in the x direction.The x distance between grid points is(6h, 4h, 2h, h, 4h, 6h).
// c
// c	Usually, there are a lot more earthquakes than seismometers, so it is
// c	generally more efficient to save travel time files from the individual
// c	stations.The program will read in an earthquake data file(header
// c	followed by a station / phase list) and then read in the appropriate
// c	travel time files
// c
// c
// c	Local Earthquakes are read in first, followed by teleseisms, and finally shots.
// c
// c	Input file attachments :
// c
// c       fn      Name   Uname   U#  Comments
// c	 1     stafile lunsta   8  Input station list
// c	 2     locdfil luneqs   9  Local Earthquake data file
// c	 3     shotfil lunsht  10  Local Shot data file(isshot / idoshot = 1)
// c	 4     telefil luntel  11  Teleseismic data file(istel / idotel = 1)
// c	 5     oldvfil luncor  12  Old coarse scale model(if ivpvs = 1)
// c	14     tabfile luntab  22  local time input files - dynamically assigned
// c	27     pbasfil lunptm  35  P base travel time file(istel / idotel = 1)
// c	28     sbasfil lunstm  36  S base travel time file(istel / idotel = 1)
// c	29     elipfil lunelp  37  Ellipticity correction  file(istel / idotel = 1)
// c
// c       Output file attachments :
// c
// c       fn      Name   Uname   U#  Comments
// c	               lout     2  sphrayderv.log(log file)
// c	 9     raystat lunray  17  Ray Statistics(if iraystat = 1)
// c	10     telrerr lunerr  18  Error summary->now in log file
// c	11     dtdsfil lundts  19  dt / ds output file
// c	12     resfile lunres  20  residual output
// c	13     hitfile lunhit  21  hit map / output
// c	15     dtdhfil lundth  23  dt / dh(hypo derivatives)
// c	16     bookfil lunbok  24  bookeeping file
// c	17     dotfile lundot  25  data out file(used if idatout = 1)
// c	18     headfil lunhed  26  output header file(used if idatout = 1)
// c	19     entfile lunmat  27  output entry points of teleseisms at base(used if nomat = 1)
// c	20     stcfile lunstc  28  Station correction derivative(used if istacor = 1)
// c	25     raypfil lunpth  33  ray path output(if iray = 1) - dynamically assigned
// c	34     sclefil lunmsc  43  Current index vector.This is used by makehyps, and is upgraded
// c	by all routines that add rows to the matrix(like addcovs, addpo, etc).
// c

#include <math.h>
#include <stdio.h>
#include "common/parseprogs.h"
#include "common/parameter.h"
#include "common/gridspec.h"
#include "common/string_process.h"
#include "common/time_process.h"
#include "common/geographic_method.h"
#include "common/environment_setting.h"
#include "common/shared_variables.h"
// #include "common/dirent.h"
#include "sphrayderv/bjdaz2.h"
#include "sphrayderv/dfind.h"
#include "sphrayderv/elpcr.h"
#include "sphrayderv/ljust.h"
#include "sphrayderv/pmin.h"

//c--mustv is the number of variables that must be assigned in the parameter
//c       file in order for this program to run correctly.The names of the
//c       variables are specified in the mvals string array below.
//c--mustf is the number of files that must be attached; see the files array below.
#define MUSTV 4
#define MUSTF 10
#define nhbyte 58 * 4

#define MAXSTRLEN 132
#define rearth 6371.0
#define degrad 0.017453292
#define hpi 3.14159265358979323846 / 2

#define _ATL_SECURE_NO_WARNINGS
#pragma warning(disable : 4206 4221 4464 4710 5045)

int kn = 0, iph = 0, kbl = 0;
int je = 0, ke = 0;

double h, dq, df;
double dsec, tarr, dpot;
double stx, sty;
float xlat, xlon;
double a, b, c, quad;
double tolmin, tolmax;
double ex, ey, ez;
double xi, yj, zk, gradt[3], dd[3], d, length, dlen, dt, fx, fy, fz, gradtm;
double sx, sy, sz, sf, sq, sr;
double xo, yo, zo;
double xold, yold, zold;
double ro, r1, r2, q1, q2, f1, f2;
double xs, ys, zs;
double sinq, cosq, tanq, ctanq, tansq, ctansq;
double dl, tc, rdevs;

float gx[nxcm], gy[nycm], gz[nzcm];
float sp[nxyzm2];
float obstime[maxobs];
float pwt[maxobs], dat[maxobs];
float rwts[maxobs], resmin[maxobs];
float tdelts[maxobs], az1s[maxobs], az2s[maxobs];
float azs[maxobs], ais[maxobs];
float tcor[maxlst][2], dtc[maxlst2];
float stlat[maxlst], stlon[maxlst];
float t[maxsta][nxyzm];
float vm[maxnbk][maxobs], vmp[maxkbl][maxobs];
float vsum[maxkbl];
float am[maxobs][4];
float avrstn[maxlst][2];
float rstnsq[maxlst][2];
float facstn[maxlst][2];
float du[nxyzm];
float slats[nxm][nym], slons[nxm][nym];

int nobstn[maxlst][2];
int ind[maxnbk][maxobs];
int mndm[maxkbl], jndx[maxmbl];
int mndex[nxyzm2], indx[nxyzm2];
int nhit[nxyzm2], mhit[nxyzm2];
int inbk[maxobs];
int jsave[maxkbl];
int istn[maxobs];
int nx, ny, nz, i, j, k, is, js, ks, ish, jsh, ksh, nseg, md, iscell, jscell, kscell, nk, nk2, nj, nj2;

char logfile[80];

char filen[80], rayfile[80];
char evid[10];
char sta[maxobs][6], stn[maxsta][6], stt[maxlst][6];
char phs[maxobs], pha[maxsta];
char mark;

// -- - the following are set so that the grid in the header can
//    be read in without overridding the vaules in the spec file

int jgridx[nxcm1], jgridy[nycm1], jgridz[nzcm1];
int wsum = 0, ncwrt = 0;
// ---- - end header stuff

int sphrayderv(char *file_parameter) {
	VERSION[10] = "2004.0923";
	char aline[MAXSTRLEN], varname[MAXSTRLEN], pval[MAXSTRLEN], parval[MAXSTRLEN],
		stafile[MAXSTRLEN], locdfil[MAXSTRLEN], shotfil[MAXSTRLEN],
		telefil[MAXSTRLEN], oldvfil[MAXSTRLEN], pbasfil[MAXSTRLEN],
		sbasfil[MAXSTRLEN], elipfil[MAXSTRLEN], raystat[MAXSTRLEN],
		telrerr[MAXSTRLEN], dtdsfil[MAXSTRLEN], resfile[MAXSTRLEN],
		hitfile[MAXSTRLEN], dtdhfil[MAXSTRLEN], bookfil[MAXSTRLEN],
		dotfile[MAXSTRLEN], headfil[MAXSTRLEN], entfile[MAXSTRLEN],
		stcfile[MAXSTRLEN], sclefil[MAXSTRLEN], specfile[MAXSTRLEN];
	int iread = 0, ivs = 1;
	double vpvs = 1.78;
	int idmean = 0, iray = 0, iraystat = 0, idatout = 1, nomat = 0;
	float resflag = 1.0;
	int ido1d = 0, ittnum = 0, ivpvs = 0, istacor = 0, idoshot = 0, idotel = 0;
	int kmin, kmax;
	int total_earthquakes = 0;
	float clat, clon;
	double xn, yn, zn;
	sscanf(file_parameter, "%s", specfile);
	specfile[MAXSTRLEN - 1] = '\0';
	FILE* fp_spc = fopen(specfile, "r");
	if (!fp_spc) {
		printf("error on opening spec-file (%s)\n", specfile);
		assert(0);
	}
	int len = 0, ierr = 0;
	char mvals[MUSTV][10] = { "nxc", "nyc", "nzc", "h" };
	for (i = 0; i < MUSTV; i++) {
		get_vars(fp_spc, mvals[i], pval, &len, &ierr);
		if (ierr == 1) {
			printf("Error trying to read variable %s", mvals[i]);
			assert(0);
		}
		if (i == 0) {
			sscanf(pval, "%d", &nxc);
		}
		else if (i == 1) {
			sscanf(pval, "%d", &nyc);
		}
		else if (i == 2) {
			sscanf(pval, "%d", &nzc);
		}
		else if (i == 3) {
			sscanf(pval, "%lf", &h);
		}
	}

	//----dimension check
	if (nxc > nxcm) {
		printf("nxc is too large.\n");
		assert(0);
	}
	if (nyc > nycm) {
		printf("nyc is too large.\n");
		assert(0);
	}
	if (nzc > nzcm) {
		printf("nzc is too large.\n");
		assert(0);
	}
	char files[MUSTF][10] = { "stafile", "locdfil", "oldvfil", "telrerr", "dtdsfil", "resfile", "hitfile", "dtdhfil", "bookfil", "sclefil" };
	char *file_list[10] = { stafile, locdfil, oldvfil, telrerr, dtdsfil, resfile, hitfile, dtdhfil, bookfil, sclefil };
	for (i = 0; i < MUSTF; i++) {
		get_vars(fp_spc, files[i], pval, &len, &ierr);
		if (ierr == 1) {
			printf("Error trying to read filename %s", files[i]);
			assert(0);
		}
		sscanf(pval, "%s", file_list[i]);
	}

	//--Optionally read in some variables
	//---Reading option
	get_vars(fp_spc, "iread ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%d", &iread);
	}
	if (iread != 0 && iread != 1) {
		iread = 0;
	}
	get_vars(fp_spc, "ivs ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%d", &ivs);
	}
	if (ivs != 0 && ivs != 1) {
		ivs = 0;
	}
	get_vars(fp_spc, "vpvs ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%lf", &vpvs);
	}
	get_vars(fp_spc, "ivpvs ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%d", &ivpvs);
	}
	if (ivpvs != 0 && ivpvs != 1) {
		ivpvs = 0;
	}

	//lzw
	float vpvsscale = 0;
	get_vars(fp_spc, "vpvsscale ", pval, &len, &ierr);
	if (ierr == 0 && ivpvs == 1) {
		sscanf(pval, "%f", &vpvsscale);
	}

	get_vars(fp_spc, "dmean ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &idmean);
	if (idmean != 0 && idmean != 1) {
		idmean = 0;
	}
	get_vars(fp_spc, "iray ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &iray);
	if (iray != 0 && iray != 1) {
		iray = 0;
	}
	get_vars(fp_spc, "iraystat ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &iraystat);
	if (iraystat != 0 && iraystat != 1) {
		iraystat = 0;
	}
	get_vars(fp_spc, "idatout ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &idatout);
	if (idatout != 0 && idatout != 1) {
		idatout = 0;
	}
	get_vars(fp_spc, "nomat ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &nomat);
	if (nomat != 0 && nomat != 1) {
		nomat = 0;
	}
	get_vars(fp_spc, "istacor ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &istacor);
	if (istacor != 0 && istacor != 1) {
		istacor = 0;
	}
	get_vars(fp_spc, "doshot ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &idoshot);
	if (idoshot != 0 && idoshot != 1) {
		idoshot = 0;
	}
	get_vars(fp_spc, "dotel ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &idotel);
	if (idotel != 0 && idotel != 1) {
		idotel = 0;
	}
	get_vars(fp_spc, "resflag ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &resflag);

	//c---Time table directory
	get_vars(fp_spc, "timedir ", pval, &len, &ierr);
	sscanf(pval, "%s", timedir);

	//c---Coordinate origin
	get_vars(fp_spc, "x0 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &x0);
	get_vars(fp_spc, "y0 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &y[0]);
	get_vars(fp_spc, "z0 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &z0);
	get_vars(fp_spc, "df ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &df);
	get_vars(fp_spc, "dq ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &dq);
	get_vars(fp_spc, "kmin ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &kmin);
	get_vars(fp_spc, "kmax ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &kmax);
	get_vars(fp_spc, "clat ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &clat);
	get_vars(fp_spc, "clon ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &clon);
	get_vars(fp_spc, "ittnum ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &ittnum);

	int ib = 0, ie = 0, lenv = 0, nvl = 0;
	rewind(fp_spc);
a11:
	aline[MAXSTRLEN - 1] = '\0';
	get_line(fp_spc, aline, &ierr);
	aline[MAXSTRLEN - 1] = '\0';
	if (ierr == 1)
		goto a12;
	if (ierr != 0)
		goto a11;
	get_field(fp_spc, aline, ib, &ie, varname, &lenv, &ierr);
	if (strncmp(varname, "igridx", lenv) != 0)
		goto a11;
	ib = ie;
	get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
	sscanf(parval, "%d", &igridx[0]);
	for (k = 1; k < nxc; k++) {
		ib = ie;
		get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
		sscanf(parval, "%d", &igridx[k]);
	}
a12:
	rewind(fp_spc);
a13:
	get_line(fp_spc, aline, &ierr);
	if (ierr == 1)
		goto a14;
	if (ierr != 0)
		goto a13;
	ib = 0;
	get_field(fp_spc, aline, ib, &ie, varname, &lenv, &ierr);
	if (strncmp(varname, "igridy", lenv) != 0)
		goto a13;
	ib = ie;
	get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
	sscanf(parval, "%d", &igridy[0]);
	for (k = 1; k < nyc; k++) {
		ib = ie;
		get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
		sscanf(parval, "%d", &igridy[k]);
	}
a14:
	rewind(fp_spc);

a15:
	get_line(fp_spc, aline, &ierr);
	if (ierr != 1) {
		if (ierr != 0)
			goto a15;
		ib = 0;
		get_field(fp_spc, aline, ib, &ie, varname, &lenv, &ierr);
		if (strncmp(varname, "igridz", lenv) != 0)
			goto a15;
		ib = ie;
		get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
		sscanf(parval, "%d", &igridz[0]);
		for (k = 1; k < nzc; k++) {
			ib = ie;
			get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
			sscanf(parval, "%d", &igridz[k]);
		}
	}

	get_vars(fp_spc, "shotfil ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", shotfil);
	get_vars(fp_spc, "telefil ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", telefil);

	get_vars(fp_spc, "pbasfil ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", pbasfil);

	get_vars(fp_spc, "sbasfil ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", sbasfil);

	get_vars(fp_spc, "elipfil ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", elipfil);

	get_vars(fp_spc, "raystat ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", raystat);

	get_vars(fp_spc, "dotfile ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", dotfile);

	get_vars(fp_spc, "headfil ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", headfil);

	get_vars(fp_spc, "entfile ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", entfile);

	get_vars(fp_spc, "stcfile ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", stcfile);

	get_vars(fp_spc, "total_earthquakes ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &total_earthquakes);

	get_vars(fp_spc, "eqkdir ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", eqkdir);
	if (eqkdir[strlen(eqkdir) - 1] != '/') {
		eqkdir[strlen(eqkdir)] = '/';
		eqkdir[strlen(eqkdir) + 1] = '\0';
	}

	//-----end of optional parameters

	nx = 1;
	ny = 1;
	nz = 1;

	for (i = 1; i < nxc; i++) {
		nx = nx + igridx[i - 1];
	}

	for (i = 1; i < nyc; i++) {
		ny = ny + igridy[i - 1];
	}

	for (i = 1; i < nzc; i++) {
		nz = nz + igridz[i - 1];
	}
	if (DEBUG_PRINT) {
		printf("  Fine grid dimension (nx, ny, nz) = %12d%12d%12d\n", nx, ny, nz);
	}
	//	c----dimension check
	if (nx > nxm) {
		printf("nx is too large.\n");
		assert(0);
	}
	if (ny > nym) {
		printf("ny is too large.\n");
		assert(0);
	}
	if (nz > nzm) {
		printf("nz is too large.\n");
		assert(0);
	}
	int nxy = nx * ny;
	int nxyz = nxy * nz;
	//c-- - Note that x0, y0, z0 are the spherical coordinates of the model origin
	//c(i.e., NW corner, top).At the start :
	//c    x0 is the geographic longitude in degress(E positive),
	//c    y0 is the geographic latitude in degrees(N positive)
	//c    z0 is depth w.r.t.to mean sea level.Positive z is down.
	//c
	//c   Note on conversion that y0 is geocentric colatitude(S positive) and z0
	//c   is still depth(so an original negative value puts the origin above
	//c   the nominal surface of the Earth(at r = 6371)
	//c   First Convert geographic latitude to geocentric colatitude
	x00 = x0;
	y00 = y[0];
	if (DEBUG_PRINT) {
		printf("  Origin:  %22.14lf %25.15lf %25.16lf       degrees/km\n", x0, y[0], z0);
	}
	y[0] *= degrad;
	y[0] = hpi - glat(y[0]);
	x0 *= degrad;
	dq *= degrad;
	df *= degrad;
	//c-- - If dq and df have not been specified, then make them so that the
	//c   interval at the surface is equal to h
	if (fabs(dq) < 0.0001)
		dq = h / rearth;
	if (fabs(df) < 0.0001)
		df = fabs(h / (rearth * sin(y[0])));

	tolmin = h * 1.E-6;
	tolmax = h * 6.0;
	if (DEBUG_PRINT) {
		printf("  Origin:  %22.14lf %25.15lf %25.16lf       radians/km\n", x0, y[0], z0);
		printf("  Radial Spacing: %24.16lf       km\n", h);
		printf("  Latitude Spacing:  %25.17E  degrees\n", dq / degrad);
		printf("  Longitude Spacing: %25.17E  degrees\n", df / degrad);
		printf("  nxc, nyc, nzc: %12d%12d%12d\n", nxc, nyc, nzc);
	}
	sprintf(logfile, "sphrayderv.log%d", ittnum);
	FILE* fp_log = fopen(logfile, "w");
	if (!fp_log) {
		printf("create %s file error.\n", logfile);
		assert(0);
	}

	fprintf(fp_log, "  \n");
	fprintf(fp_log,
		" *************************************************************** \n");
	fprintf(fp_log, "          Parameters Set For This Run of Sphrayderv \n");
	fprintf(fp_log, "  \n");
	fprintf(fp_log, " Sphrayderv VERSION: %s\n", VERSION);
	fprintf(fp_log, "  \n");
	fprintf(fp_log, " Current parameter specification file: %-40s\n", specfile);
	fprintf(fp_log, "  \n");
	fprintf(fp_log, "  Iteration counter: %18d\n", ittnum);
	fprintf(fp_log, "  \n");
	fprintf(fp_log, "  Longitude origin (x0) : %12.5lf    \n", x00);
	fprintf(fp_log, "  Latitude origin  (y0) : %12.5lf    \n", y00);
	fprintf(fp_log, "  Depth origin     (z0) : %21.16lf     \n", z0);
	{
		char tmp[MAXSTRLEN];
		dtoa(tmp, h, 18);
		fprintf(fp_log, "  Fine Radial Spacing   :   %s       km\n", tmp);
	}
	fprintf(fp_log, "  Fine Latitude Spacing :   %19.17E  degrees\n", dq / degrad);
	fprintf(fp_log, "  Fine Longitude Spacing:   %19.17E  degrees\n", df / degrad);
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of X coarse grid nodes: %12d\n", nxc);
	fprintf(fp_log, "  X coarse grid node spacing: \n");
	for (i = 0; i < nxc - 1; i++) {
		fprintf(fp_log, "% 4d", igridx[i]);
		if (i % 10 == 9) {
			fprintf(fp_log, "\n");
		}
	}
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of Y coarse grid nodes: %12d\n", nyc);
	fprintf(fp_log, "  Y coarse grid node spacing: \n");
	for (i = 0; i < nyc - 1; i++) {
		fprintf(fp_log, "% 4d", igridy[i]);
		if (i % 10 == 9) {
			fprintf(fp_log, "\n");
		}
	}
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of Z coarse grid nodes: %12d\n", nzc);
	fprintf(fp_log, "  Z coarse grid node spacing: \n");
	for (i = 0; i < nzc - 1; i++) {
		fprintf(fp_log, "% 4d", igridz[i]);
		if (i % 10 == 9) {
			fprintf(fp_log, "\n");
		}
	}
	fprintf(fp_log, "\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of X fine grid nodes: %12d\n", nx);
	fprintf(fp_log, "  Number of Y fine grid nodes: %12d\n", ny);
	fprintf(fp_log, "  Number of Z fine grid nodes: %12d\n", nz);
	fprintf(fp_log, " \n");

	fprintf(fp_log, " Travel Time Table Directory: %-40s\n", timedir);
	fprintf(fp_log, " \n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Input file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Input station list: %-53s\n", stafile);
	fprintf(fp_log, " Local Earthquake data file: %-31s\n", locdfil);
	if (idoshot)
		fprintf(fp_log, " Local Shot data file: %s\n", shotfil);

	if (idotel) {
		fprintf(fp_log, " Teleseismic data file: %s\n", telefil);
		fprintf(fp_log, " P base travel time file: %s\n", pbasfil);
		fprintf(fp_log, " S base travel time file: %s\n", sbasfil);
		fprintf(fp_log, " Ellipticity correction file: %s\n", elipfil);
	}
	if (ivpvs)
		fprintf(fp_log, " Old coarse scale model: %s\n", oldvfil);
	if (istacor)
		fprintf(fp_log, " Station correction dervs: %s\n", stcfile);

	fprintf(fp_log, " \n");
	fprintf(fp_log, " Output file attachments:\n");
	fprintf(fp_log, " \n");
	if (iraystat)
		fprintf(fp_log, " Ray Statistics file: %s\n", raystat);
	fprintf(fp_log, " Error Summary file: %s\n", telrerr);
	fprintf(fp_log, " dt/ds file: %s\n", dtdsfil);
	fprintf(fp_log, " Residuals file: %s\n", resfile);
	fprintf(fp_log, " Hit file: %s\n", hitfile);
	fprintf(fp_log, " Hypo derivatives file: %s\n", dtdhfil);
	fprintf(fp_log, " Book keeping file: %s\n", bookfil);
	if (idatout) {
		fprintf(fp_log, " Data Output file: %s\n", dotfile);
		fprintf(fp_log, " Header Output file: %s\n", headfil);
	}
	if (nomat) {
		fprintf(fp_log, " Tele Entry Point file: %s\n", entfile);
		fprintf(fp_log, " Current index vector: %s\n", sclefil);
	}
	fprintf(fp_log, " \n");

	int nxyc = nxc * nyc;
	int nxyzc = nxyc * nzc;
	int nxyzc2 = nxyzc * 2;

	double xmax = x0 + (nx - 1) * df;
	double ymax = y[0] + (ny - 1) * dq;
	double zmax = z0 + (nz - 1) * h;

	if (DEBUG_PRINT) {
		printf("  xmax = %21.14lf       degrees\n", xmax / degrad);
		printf("  ymax = %21.14lf       degrees\n", ymax / degrad);
		printf("  zmax = %13.6lf      km\n", zmax);
	}

	fprintf(fp_log, "  Total Number of fine grid nodes:   %12d\n", nxyz);
	fprintf(fp_log, "  Total Number of coarse grid nodes: %12d\n", nxyzc);
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  X Max:    %18.14lf       degrees\n", xmax / degrad);
	fprintf(fp_log, "  Y Max:    %18.14lf       degrees(colatitude)\n", ymax / degrad);
	fprintf(fp_log, "  Y Max:    %18.14lf       degrees(latitude)\n", 90. - ymax / degrad);
	fprintf(fp_log, "  Z Max: %12.6lf      km\n", zmax);
	fprintf(fp_log, " \n");

	char doo[2][10] = { " not", "" };
	fprintf(fp_log, "  Sphrayderv will%s do teleseisms (idotel = %d)\n", doo[idotel], idotel);
	fprintf(fp_log, "  Sphrayderv will%s do shots (idoshot = %d)\n", doo[idoshot], idoshot);
	fprintf(fp_log, "  Sphrayderv will%s demean shot residuals (idmean=%d)\n", doo[idmean], idmean);
	fprintf(fp_log, "  Station corrections are%s variable (istacor = %d)\n", doo[istacor], istacor);
	if (ivpvs) {
		fprintf(fp_log, "  Vp/Vs is solved for (ivpvs = 1)\n");
		//clzw
		fprintf(fp_log, "  Vp/Vs is scaling by vpvsscale=%f", vpvsscale);
	}
	else {
		fprintf(fp_log, "  Shear slowness is solved for (ivpvs = 0)\n");
	}
	strcpy(doo[0], " does not");
	strcpy(doo[1], "");
	fprintf(fp_log, "  Sphrayderv%s output raypaths (iray = %d)\n", doo[iray], iray);
	fprintf(fp_log, "  Sphrayderv%s output ray stats (iraystat = %d)\n", doo[iraystat], iraystat);

	fprintf(fp_log, "  Sphrayderv%s creates data file (idatout = %d)\n", doo[idatout], idatout);
	if (idatout)
		fprintf(fp_log, "  Residuals are flagged if more than %12.7f\n", resflag);
	fprintf(fp_log, "  Sphrayderv %s make matrices (nomat = %d)\n", doo[!nomat], nomat);
	strcpy(doo[0], "3D");
	strcpy(doo[1], "1D");
	fprintf(fp_log, "  Sphrayderv will do a %s inversion (do1d = %d)\n", doo[ido1d], ido1d);
	fprintf(fp_log, " \n");

	//c---- - m1 and m2 are used in bookkeeping to notify the makehyps program of the presence of a "non-earthquake"
	//c     which requires an adjustment only to origin time(in the case of m1) or nothing at all(m2).Often
	//c     we choose m1 for shots since there may be some uncertainty in the shot time or local wavespeeds, but
	//c     we can optionally make this m2 if we wish to match absolute times.
	//c
	//c     Note that at present this is hardwired(shots wired to m2)
	int m1 = -1, m2 = -2;
	FILE* fp_sta = fopen(stafile, "r");
	if (!fp_sta) {
		printf("Error on opening file: %s\n", stafile);
		assert(0);
	}
	FILE* fp_eqs = fopen(locdfil, "r");
	if (!fp_eqs) {
		printf("Error on opening file: %s\n", locdfil);
		assert(0);
	}
	FILE* fp_ray = NULL;
	if (iraystat) {
		fp_ray = fopen(raystat, "w");
		if (!fp_ray) {
			printf("create %s file error.\n", raystat);
			assert(0);
		}
	}
	FILE* fp_err = fopen(telrerr, "w");
	if (!fp_err) {
		printf("create %s file error.\n", telrerr);
		assert(0);
	}
	FILE* fp_dts = fopen(dtdsfil, "wb");
	if (!fp_dts) {
		printf("create %s file error.\n", dtdsfil);
		assert(0);
	}
	FILE* fp_dth = fopen(dtdhfil, "wb");
	if (!fp_dth) {
		printf("create %s file error.\n", dtdhfil);
		assert(0);
	}
	FILE* fp_res = fopen(resfile, "wb");
	if (!fp_res) {
		printf("create %s file error.\n", resfile);
		assert(0);
	}
	FILE* fp_hit = fopen(hitfile, "wb");
	if (!fp_hit) {
		printf("create %s file error.\n", hitfile);
		assert(0);
	}
	FILE* fp_bok = fopen(bookfil, "wb");
	if (!fp_bok) {
		printf("create %s file error.\n", bookfil);
		assert(0);
	}
	FILE* fp_msc = fopen(sclefil, "wb");
	if (!fp_msc) {
		printf("create %s file error.\n", sclefil);
		assert(0);
	}

	FILE* fp_stc = NULL;
	FILE* fp_mat = NULL;
	FILE *fp_dot = NULL, *fp_hed = NULL;
	/*
	tabfile(local time input files - dynamically assigned)		.. / data / small / TTimes00
	locdfil(Local Earthquake data file)							.. / data / small / local.data_re01
	stafile(Input station list)									.. / data / small / runs_files / stationloc_out.txt

	oldvfil(Old coarse scale model(if ivpvs = 1))				.. / data / small / TW_m30_wL.mod

	shotfil(Local Shot data file(isshot / idoshot = 1))			.. / data / small / runs_files / arrivals / all_shot.data_04
	telefil(Teleseismic data file(istel / idotel = 1))			.. / data / small / runs_files / arrivals / tele_01_09.dat
	pbasfil(P base travel time file(istel / idotel = 1))		.. / data / small / runs_files / parameter / Ptimes.base293
	sbasfil(S base travel time file(istel / idotel = 1))		.. / data / small / runs_files / parameter / Stimes.base293
	elipfil(Ellipticity correction  file(istel / idotel = 1))	.. / data / small / runs_files / parameter / elpcor
	*/

	if (idatout) {
		fp_dot = fopen(dotfile, "w");
		if (!fp_dot) {
			printf("create %s file error.\n", dotfile);
			assert(0);
		}
		fp_hed = fopen(headfil, "w");
		if (!fp_hed) {
			printf("create %s file error.\n", headfil);
			assert(0);
		}
	}

	if (istacor) {
		fp_stc = fopen(stcfile, "wb");
		if (!fp_stc) {
			printf("create %s file error.\n", stcfile);
			assert(0);
		}
	}

	if (nomat) {
		fp_mat = fopen(entfile, "wb");
		if (!fp_mat) {
			printf("create %s file error.\n", entfile);
			assert(0);
		}
	}

	gx[0] = x0;
	gy[0] = y[0];
	gz[0] = z0;
	for (i = 1; i < nxc; i++) {
		gx[i] = gx[i - 1] + df * igridx[i - 1];
	}
	for (i = 1; i < nyc; i++) {
		gy[i] = gy[i - 1] + dq * igridy[i - 1];
	}
	for (i = 1; i < nzc; i++) {
		gz[i] = gz[i - 1] + h * igridz[i - 1];
	}

	int nxyz2 = nxyz * 2;
	int maxvar = nxyz;
	int maxvar2 = nxyz2;
	int maxvarc = nxyzc;
	int maxvarc2 = nxyzc2;
	if (ido1d) {
		maxvar = nz;
		maxvar2 = nz * 2;
		maxvarc = nzc;
		maxvarc2 = nzc * 2;
	}

	int nstr;
	char str_inp[MAXSTRLEN];
	double stz[maxlst];

	for (nstr = 0; fgets(str_inp, sizeof(str_inp), fp_sta); nstr++) {
		char str_inp_trimed[MAXSTRLEN];
		trimwhitespace(str_inp_trimed, strlen(str_inp), str_inp);
		if (str_inp_trimed[0] == '\n')
			break;
		if (strlen(str_inp_trimed) > MAXSTRLEN) {
			printf("str_inp is too long: %s\n", str_inp_trimed);
			assert(0);
		}
		if (nstr >= maxlst) {
			printf(" Error: too many stations in list.. aborting\n");
			fprintf(fp_log, " Error: too many stations in list.. aborting\n");
			assert(0);
		}
		//
		// ---- sty and stx are in meters from the origin when using cartesian coordinates,
		//    but they are not used in this version.  slat and slon are lats
		//    and lons in decimal degrees.  stx and sty save the lats and lons relative
		//    to the origin.  ielev is elevation in meters above m.s.l.
		//

		int ielev;
		double slat, slon;
		sscanf(str_inp, "%lf %lf %d %s %lf %lf %f %f", &sty, &stx, &ielev, stt[nstr], &slat, &slon, &tcor[nstr][0], &tcor[nstr][1]);
		// ---temp lines to reassign slat and slon as the reference locations
		stx = slon;
		sty = slat;
		// --end temp lines
		stlon[nstr] = stx * degrad;
		stlat[nstr] = hpi - glat(sty * degrad);
		stz[nstr] = -ielev / 1000.;
		//    stz(nstr) = -elev
		// ---left justify the station name (no leading blanks)
		ljust(stt[nstr]);
	}
	printf(" %d stations read in.\n", nstr);
	fclose(fp_sta);

	// -- - header on station correction array
	int nstr2 = 2 * nstr;
	if (istacor) {
		fprintf(fp_stc, " %d", nstr2);
	}

	for (i = 0; i < nstr; i++) {
		for (j = 0; j < 2; j++) {
			avrstn[i][j] = 0;
			rstnsq[i][j] = 0;
			facstn[i][j] = 0;
			nobstn[i][j] = 0;
		}
	}

	// ---- - set the data file to the local earthquake data file
	FILE* fp_din = fp_eqs;
	int mbl = 0;
	for (int i1 = 0; i1 < nxyz2; i1++) {
		nhit[i1] = 0;
	}

	int nev = 0, nshot = 0, ntel = 0, ntread = 0, iover = 0;
	int knobs = 0, isshot = 0, istel = 0;
	float facs = 0, facsev = 0, avres = 0, avresev = 0;
	float rms = 0, rsq = 0, wtmod = 0;
	float xm = 0, ym = 0, zm = 0;

	// ****************** Start Loop over events ******************
	//    Read in event header
a3:
	fgets(str_inp, sizeof(str_inp), fp_din);
	if (feof(fp_din)) {
		goto a60;
	}
	if (strlen(str_inp) > sizeof(str_inp)) {
		printf("str_inp is too long: %s\n", str_inp);
		assert(0);
	}
	int kyr, kday, khr, kmn;
	double esec;
	float dep;

	sscanf(str_inp, "%d %d %d %d %lf %f %f %f %s\n", &kyr, &kday, &khr,
		&kmn, &esec, &xlat, &xlon, &dep, evid);

	//	c---convert to epochal time
	dsec = esec;
	htoe(kyr, kday, khr, kmn, dsec, &dpot);

	if (isshot == 0) {
		if (istel == 0) {
			nev++;
		}
		else {
			ntel++;
		}
	}
	else {
		nshot++;
	}
	int ies = 0, jes = 0, kes = 0;
	double xis = 0, yjs = 0, zks = 0;
	if (istel == 0) {
		ex = xlon * degrad;
		ey = hpi - glat(xlat * degrad);
		ez = dep;

		// -- - Grid point 1 for the earthquake or shot
		ies = (int)((ex - x0) / df);
		jes = (int)((ey - y[0]) / dq);
		kes = (int)((ez - z0) / h);

		// -- - make sure event is within the model
		if ((ex < x0 || ies >= nx - 1)
			|| (ey < y[0] || jes >= ny - 1)
			|| (ez < z0 || kes >= nz - 1)) {
			printf(" Error:  Earthquake out of bounds; skipping..\n");
			fprintf(fp_err, " Error:  This event is out of bounds and skipped : \n");
			fprintf(fp_err, " %4d %3d %2d %2d %8.4lf %9.5f %10.5f %8.4lf %12s\n", kyr, kday, khr, kmn, esec, ex, ey, ez, evid);
			fprintf(fp_err, "\n");
			while (aline[0] != '\0') {
				if(fgets(str_inp, sizeof(str_inp), fp_din)== NULL) {
					break;
				}
				len = (int)strlen(str_inp);
				if (len > MAXSTRLEN) {
					printf("input length is too large. len=%d str_inp=%s\n", len,
						str_inp);
					assert(0);
				}
				strcpy(str_inp, trim(str_inp));
				if (str_inp[0] == '\n' || str_inp[0] == '\0') {
					break;
				}
			}
			goto a3;
		}

		// ----coordinates of point 1
		xis = df * ies + x0;
		yjs = dq * jes + y[0];
		zks = h * kes + z0;
	}

	// ----Read in the times and stations
	avresev = 0;
	facsev = 0;
	int nsta = 0;
	int isgood[maxobs];
	memset(isgood, 0, sizeof(isgood));
	for (nsta = 0; fgets(str_inp, sizeof(str_inp), fp_din); nsta++) {
		if (nsta >= maxobs) {
			printf("Error:  too many observations! nsta=%d maxobs=%d\n", nsta, maxobs);
			fprintf(fp_log, "Error:  too many observations! nsta=%d maxobs=%d\n", nsta, maxobs);
			assert(0);
		}
		char str_inp_trimed[100];
		trimwhitespace(str_inp_trimed, strlen(str_inp), str_inp);

		len = (int)strlen(str_inp_trimed);
		if (str_inp_trimed[0] == '\0') {
			break;
		}
		if (len >= MAXSTRLEN) {
			printf("input length is too large. len=%d str_inp=%s\n", len,
				str_inp_trimed);
			assert(0);
		}

		char str_tmp[MAXSTRLEN];
		int iyr, jday, ihr, imn;
		double sec = 0;
		sscanf(str_inp_trimed, "%s %d %d %d %d %lf %99[^\n]\n", sta[nsta], &iyr,
			&jday, &ihr, &imn, &sec, str_tmp);
		if (str_tmp[0] == '*') {
			str_tmp[0] = ' ';
			//trim(str_tmp);
		}
		char str_t[10];
		sscanf(str_tmp, " %c %f", str_t, &rwts[nsta]);
		phs[nsta] = str_t[0];
		// --- left justify station name
		ljust(sta[nsta]);
		// ---convert from "0 1" to "P S" if necessary
		if (phs[nsta] == '0') {
			phs[nsta] = 'P';
		}
		else if (phs[nsta] == '1') {
			phs[nsta] = 'S';
		}
		d_blank(filen, &len);

		pwt[nsta] = 1. / (rwts[nsta] * rwts[nsta]);
		dsec = sec;

		htoe(iyr, jday, ihr, imn, dsec, &tarr);
		obstime[nsta] = tarr - dpot;
		if (nsta > maxobs) {
			printf(" Error: too many observations!\n");
			fprintf(fp_log, " Error: too many observations!\n");
			assert(0);
		}
	}
	if (DEBUG_PRINT) {
		if (isshot == 0) {
			if (istel == 0) {
				printf(" %d times read in for local event %d\n", nsta, nev);
				printf(" Working on local event %s\n", evid);
			}
			else {
				printf(" %d times read in for teleseismic event %d\n", nsta, ntel);
				printf(" Working on tele event %s\n", evid);
			}
		}
		else {
			printf(" %d times read in for shot %d\n", nsta, nshot);
			printf(" Working on shot %s\n", evid);
		}
	}

	// *****Start Loop over Phases*****
	for (i = 0; i < nsta; i++) {

		// ------clear derivative array
		memset(du, 0, maxvar * sizeof(du[0]));

		inbk[i] = 0;
		rdevs = 0;
		isgood[i] = 1;

		// Find this station in the station list
		for (kn = 0; kn < nstr; kn++) {
			if (strcmp(sta[i], stt[kn]) == 0) {
				break;
			}
		}
		if (kn == nstr) {
			printf(" Warning: Station not in list = %s ...skipping this phase.\n", sta[i]);
			fprintf(fp_log, " Warning: Station not in list = %s ...skipping this phase.\n", sta[i]);
			isgood[i] = 0;
			goto a200;
			// continue; // for(i=0;i<nsta;i++) {...}
		}
		xs = stlon[kn];
		ys = stlat[kn];
		zs = stz[kn];

		iph = 0;
		if (phs[i] == 'S')
			iph = 1;
		int iphm1 = iph;
		tc = tcor[kn][iph];
		istn[i] = kn + nstr * iphm1;

		//--- This is the nearest grid point to the Station
		is = (int)(round((xs - x0) / df));
		js = (int)(round((ys - y[0]) / dq));
		ks = (int)(round((zs - z0) / h));

		//--- This is the cell containing the Station
		iscell = (int)((xs - x0) / df);
		jscell = (int)((ys - y[0]) / dq);
		kscell = (int)((zs - z0) / h);
		// printf("1083 iscell=%d\tjsecll=%d\tkscell=%d\n", iscell, jscell, kscell);
		// stop

		//---See if the travel time tables for this station has been read in.  If not, read them in.
		if (ntread > 0) {
			for (j = 0; j < ntread; j++) {
				if (strcmp(sta[i], stn[j]) == 0) {
					if (ivs == 0 || phs[i] == pha[j]) {
						goto a7;
					}
				}
			}
		}

		int iuse = ntread;
		ntread++;
		if (ntread > maxsta) {
			if (DEBUG_PRINT)
				printf(" Limit reached...looking for one we do not need now.\n");
			assert(0);
			// ---Loop over the station-phases that have been read in.
			// If it happens that one of these is not used by the current event, read over it.
			// If there is no more space available, choose one at random to read over.

			for (j = 0; j < ntread - 1; j++) {
				int flag = 0;
				for (k = 0; k < nsta; k++) {
					if (strcmp(sta[k], stn[j]) == 0) {
						if (ivs == 0 || phs[k] == pha[j]) {
							flag = 1;
							break;
						}
					}
				}
				if (flag)
					continue;
				goto a21;
			}
			iover++;
			if (iover > maxsta)
				iover = 1;
			iuse = iover;
			j = iuse;
		a21:
			ntread--;
			iuse = j;
			if (DEBUG_PRINT)
				printf("Overwriting %s  %c with  %s %c\n", stn[j], pha[j], sta[i], phs[i]);
		}

		strcpy(stn[iuse], sta[i]);
		pha[iuse] = phs[i];
		if (phs[i] == 'P' || ivs == 0) {
			sprintf(filen, "%s/%s.ptimes", timedir, sta[i]);
		}
		else {
			sprintf(filen, "%s/%s.stimes", timedir, sta[i]);
		}

		d_blank(filen, &len);
		printf("Reading in %s\n", filen);

		FILE* fp_tab = fopen(filen, "rb");
		if (!fp_tab) {
			printf("Error on opening file: %s\n", filen);
			assert(0);
		}

		char hdr[nhbyte];
		char head[5], type[5], syst[5];
		char quant[5];
		char flatten[5];
		char hcomm[125];

		fread(hdr, sizeof(hdr[0]), nhbyte, fp_tab);
		char* offset = hdr;
		sscanf(offset, "%4s", head);
		offset += strlen(head);
		sscanf(offset, "%4s", type);
		offset += strlen(type);
		sscanf(offset, "%4s", syst);
		offset += strlen(syst);
		sscanf(offset, "%4s", quant);
		offset += strlen(quant);
		sscanf(offset, "%4s", flatten);
		offset += strlen(flatten);
		sscanf(offset, "%124s", hcomm);

		// ---see if this the time file has a header on it.
		if (strcmp(head, "HEAD") == 0) {
			if (DEBUG_PRINT) {
				printf(" File has a header...\n");
			}
			if (strcmp(type, "FINE") != 0) {
				if (DEBUG_PRINT)
					printf("WARNING: input mesh does not appear to be FINE: %s\n", type);
			}
			if (iread == 0) {
				fread(t[iuse], sizeof(t[iuse][0]), nxyz, fp_tab);
			}
		}
		else {
			rewind(fp_tab);
			if (DEBUG_PRINT) {
				printf(" File has no header...\n");
			}
			if (iread == 0) {
				fread(t[iuse], sizeof(t[iuse][0]), nxyz, fp_tab);
			}
		}
		fclose(fp_tab);
		if (DEBUG_PRINT)
			printf("...Done.\n");
		j = iuse;

	a7:
		// ----test to see if all data can be read in correctly
		if (iread == 1) {
			goto a3;
		}

		// ----(rsx, rsy) is a unit vector from earthquake to station used in raypath stats only
		double rsx = xs - ex;
		double rsy = ys - ey;

		// ----if this is a teleseism, find the entry point at the base of the model
		double ttel = 0, ebx = 0, eby = 0, ebz = 0;
		if (istel) {
			// --- calculate the ellipticity correction for this station - event pair(need station at surface!)
			double rlon = xlon * degrad;
			double rlat = hpi - glat(xlat * degrad);
			ez = dep;
			float slat = stlat[kn];
			float slon = stlon[kn];
			float delt = 0, az1 = 0, az2 = 0;
			bjdaz2(slat, slon, rlat, rlon, &delt, &az1, &az2, 0);
			double tdelt = delt / degrad;
			tdelts[i] = tdelt;
			az1s[i] = az1 / degrad;
			az2s[i] = az2 / degrad;
			// ------elpcr gets the proper ellipticity correction for this ray
			double elat = hpi - rlat;
			double azdd = az2 / degrad;
			double ezr4 = ez;
			double elpc = 0;
			elpcr(elat, tdelt, azdd, ezr4, &elpc, iphm1, 1);
			double time = 0, tbase = 0;
			int inbound = 0;
			// entry(rlat, rlon, ezr4, &ebx, &eby, &ebz, iphm1, j, elpc, &time, &tbase, &inbound);

			if (inbound == 0) {
				printf(" Error: this station has an out of bounds start point: %s\n", sta[i]);
				printf(" Skipping this station\n");
				//write(lunerr, 102) kyr, kday, khr, kmn, esec, ex, ey, ez, evid
				printf(" Error:  this station has an out of bounds start point: %s delt=%lf\n", sta[i], tdelt);
				printf("\n");
				isgood[i] = 0;
				goto a200;
				// continue; // for(i=0;i<nsta;i++){...
			}
			if (nomat) {
				double xent = ebx / degrad;
				double yent = glatinv(hpi - eby) / degrad;
				fprintf(fp_mat, "%lf %lf %lf\n", xent, yent, ebz);
				// --- extra output for testing
			}
			ttel = tbase;

			// --- Grid point 1 for the entry point
			ies = (int)((ebx - x0) / df);
			jes = (int)((eby - y[0]) / dq);
			kes = (int)((ebz - z0) / h);
			if (((ebx - x0) < 0.) || (ies >= nx) ||
				((eby - y[0]) < 0.) || (jes >= ny) ||
				((ebz - z0) < 0.) || (kes >= nz)) {
				printf(" Error:  this station has an out of bounds entry point: %s\n", sta[i]);
				printf(" Skipping this station.\n");
				fprintf(fp_err, " %d %d %d %d %lf %lf %lf %lf %s\n", kyr, kday, khr, kmn, esec, ex, ey, ez, evid);
				printf(" Error:  this station has an out of bounds entry point: %s\n", sta[i]);
				printf("\n");
				isgood[i] = 0;
				goto a200;
				// continue; // for(i=0;i<nsta;i++){...
			}

			// ----coordinates of point 1
			xis = df * ies + x0;
			yjs = dq * jes + y[0];
			zks = h * kes + z0;

			rsx = xs - ebx;
			rsy = ys - eby;
		}
		double smag = sqrt(rsx * rsx + rsy * rsy);
		rsx /= smag;
		rsy /= smag;

		// ----- Now Trace ray from the earthquake to the station
		nseg = 0;
		int iseg = 0;
		int isegm1 = 0;
		int isegm2 = 0;

		// ----start at the earthquake
		double xx = 0, yy = 0, zz = 0;
		if (istel == 0) {
			xx = ex;
			yy = ey;
			zz = ez;
		}
		else {
			xx = ebx;
			yy = eby;
			zz = ebz;
		}
		// --- (ie, je, ke) is point 1 for the current cell
		ie = ies;
		je = jes;
		ke = kes;

		// --- coordinates of point 1
		xi = xis;
		yj = yjs;
		zk = zks;

		// --- Calculate the traveltime at the earthquake by trilinear interpolation

		nk = nxy * ke;
		nj = nx * je;
		nk2 = nxy * (ke + 1);
		nj2 = nx * (je + 1);

		//c Point 1 = SW TOP = nk + nj + ie
		//c Point 2 = SE TOP = nk + nj + ie + 1
		//c Point 3 = NW TOP = nk + nj2 + ie
		//c Point 4 = NE TOP = nk + nj2 + ie + 1
		//c Point 5 = SW BOT = nk2 + nj + ie
		//c Point 6 = SE BOT = nk2 + nj + ie + 1
		//c Point 7 = NW BOT = nk2 + nj2 + ie
		//c Point 8 = NE BOT = nk2 + nj2 + ie + 1

		int ipt[8];
		ipt[0] = nk + nj + ie;
		ipt[1] = nk + nj + ie + 1;
		ipt[2] = nk + nj2 + ie;
		ipt[3] = nk + nj2 + ie + 1;
		ipt[4] = nk2 + nj + ie;
		ipt[5] = nk2 + nj + ie + 1;
		ipt[6] = nk2 + nj2 + ie;
		ipt[7] = nk2 + nj2 + ie + 1;

		// Don't use this ray if it's in a region of the model where times are not
		// alculated(see parameter "maxoff" in punch.c)

		for (int iii = 0; iii < 8; iii++) {
			if (t[j][ipt[iii]] > 1.E9) {
				printf(" Ray out of bounds...skipping this observation\n");
				fprintf(fp_err, " Ray out of bounds...skipping this observation\n");
				isgood[i] = 0;
			}
		}
		if (isgood[i] == 0) {
			goto a200;
			// continue; // for(i=0;i<nsta;i++){...
		}
		// printf("1345 xx=%f xi=%f df=%f\n", xx, xi, df);
		fx = (xx - xi) / df;
		fy = (yy - yj) / dq;
		fz = (zz - zk) / h;
		// printf("1345 fx=%lf fy=%lf fz=%lf\n", fx, fy, fz);

		double ds[8];
		ds[0] = (1. - fz)*(1. - fy)*(1. - fx);
		ds[1] = (1. - fz)*(1. - fy)*fx;
		ds[2] = (1. - fz)*fy*(1. - fx);
		ds[3] = (1. - fz) * fy * fx;
		ds[4] = fz * (1. - fy) * (1. - fx);
		ds[5] = fz * (1. - fy) * fx;
		ds[6] = fz * fy * (1. - fx);
		ds[7] = fz * fy * fx;

		dt = 0;
		for (int iii = 0; iii < 8; iii++) {
			dt += ds[iii] * t[j][ipt[iii]];
		}

		// --- For local events, ttel = 0. For teles, ttel is the time to the base of the model.
		dt += ttel;

		// --- apply station correction
		dt += tc;
		dt = obstime[i] - dt;
		resmin[i] = dt;

		//     This could be skipped for teles, but is innocuous
		double dsdx[8];
		dsdx[0] = -(1. - fz) * (1. - fy);
		dsdx[1] = (1. - fz) * (1. - fy);
		dsdx[2] = -(1. - fz) * fy;
		dsdx[3] = (1. - fz) * fy;
		dsdx[4] = -fz * (1. - fy);
		dsdx[5] = fz * (1. - fy);
		dsdx[6] = -fz * fy;
		dsdx[7] = fz * fy;

		double dsdy[8];
		dsdy[0] = -(1. - fz) * (1. - fx);
		dsdy[1] = -(1. - fz) * fx;
		dsdy[2] = (1. - fz) * (1. - fx);
		dsdy[3] = (1. - fz) * fx;
		dsdy[4] = -fz * (1. - fx);
		dsdy[5] = -fz * fx;
		dsdy[6] = fz * (1. - fx);
		dsdy[7] = fz * fx;

		double dsdz[8];
		dsdz[0] = -(1. - fy) * (1. - fx);
		dsdz[1] = -(1. - fy) * fx;
		dsdz[2] = -fy * (1. - fx);
		dsdz[3] = -fy * fx;
		dsdz[4] = (1. - fy) * (1. - fx);
		dsdz[5] = (1. - fy) * fx;
		dsdz[6] = fy * (1. - fx);
		dsdz[7] = fy * fx;
		double dtdx = 0, dtdy = 0, dtdz = 0;
		for (int iii = 0; iii < 8; iii++) {
			dtdx += dsdx[iii] * t[j][ipt[iii]];
			dtdy += dsdy[iii] * t[j][ipt[iii]];
			dtdz += dsdz[iii] * t[j][ipt[iii]];
		}
		dtdx /= df;
		dtdy /= dq;
		dtdz /= h;

		//----- wt = 1 / variance
		//----- wet = 1 / standard deviation
		double wt = pwt[i];
		double wet = sqrt(wt);

		double res1 = dt * wt;
		double res2 = res1 * dt;
		avres += res1;
		avresev += res1;
		rsq += res2;
		rms += res1 * res1;
		facs += wt;
		facsev += wt;
		wtmod += wt * wt;
		knobs++;

		// ----collect stats for station
		avrstn[kn][iph] += res1;
		rstnsq[kn][iph] += res2;
		facstn[kn][iph] += wt;
		nobstn[kn][iph]++;

		// ***** Propagate ray from earthquake to station ******

		//----location of cell
		double dxs = 0, dys = 0, dzs = 0;
		double rs = 0, dqs = 0, dfs = 0;
	a100:
		xi = df * ie + x0;
		yj = dq * je + y[0];
		zk = h * ke + z0;

		//     Test to see if ray has reached the station.If so, accumulate normal equation
		//     stuff and then go do next ray

		dxs = xx - xs;
		dys = yy - ys;
		dzs = zz - zs;
		rs = rearth - zz;
		dqs = rs * dys;
		dfs = rs * sin(yy) * dxs;
		dlen = sqrt(dfs * dfs + dqs * dqs + dzs * dzs);

		if ((ie == iscell && je == jscell && ke == kscell) || (dlen < (h / 1000.))) {

			// ----add on derivative for the station cell segment
			xm = xx - dxs * 0.5;
			ym = yy - dys * 0.5;
			zm = zz - dzs * 0.5;
			int i1 = -1, j1 = -1, k1 = -1;

			//--- find xm in coarse grid
			for (int iii = 1; iii < nxc; iii++) {
				if (xm < gx[iii]) {
					i1 = iii;
					break;
				}
			}
			if (i1 == -1) {
				printf(" Error xm = %lf is out of bounds!\n", xm);
				fprintf(fp_err, " Error xm = %lf is out of bounds!\n", xm);
				fprintf(fp_log, " Error xm = %lf is out of bounds!\n", xm);
				assert(0);
			}
			for (int iii = 1; iii < nyc; iii++) {
				if (ym < gy[iii]) {
					j1 = iii;
					break;
				}
			}
			if (j1 == -1) {
				printf(" Error ym = %lf is out of bounds!\n", ym);
				fprintf(fp_err, " Error ym = %lf is out of bounds!\n", ym);
				fprintf(fp_log, " Error ym = %lf is out of bounds!\n", ym);
				assert(0);
			}
			for (int iii = 1; iii < nzc; iii++) {
				if (zm < gz[iii]) {
					k1 = iii;
					break;
				}
			}
			if (k1 == -1) {
				printf(" Error zm = %lf is out of bounds!\n", zm);
				fprintf(fp_err, " Error zm = %lf is out of bounds!\n", zm);
				fprintf(fp_log, " Error zm = %lf is out of bounds!\n", zm);
				assert(0);
			}
			float hx = gx[i1] - gx[i1 - 1];
			float hy = gy[j1] - gy[j1 - 1];
			float hz = gz[k1] - gz[k1 - 1];
			float xci = gx[i1 - 1];
			float ycj = gy[j1 - 1];
			float zck = gz[k1 - 1];

			//--- derivatives
			fx = (xm - xci) / hx;
			fy = (ym - ycj) / hy;
			fz = (zm - zck) / hz;

			ds[0] = (1. - fz) * (1. - fy) * (1. - fx);
			ds[1] = (1. - fz) * (1. - fy) * fx;
			ds[2] = (1. - fz) * fy * (1. - fx);
			ds[3] = (1. - fz) * fy * fx;
			ds[4] = fz * (1. - fy) * (1. - fx);
			ds[5] = fz * (1. - fy) * fx;
			ds[6] = fz * fy * (1. - fx);
			ds[7] = fz * fy * fx;

			nk = nxyc * (k1 - 1);
			nj = nxc * (j1 - 1);
			nk2 = nxyc * k1;
			nj2 = nxc * j1;

			int i1m1 = i1;

			ipt[0] = nk + nj + i1m1;
			ipt[1] = nk + nj + i1m1 + 1;
			ipt[2] = nk + nj2 + i1m1;
			ipt[3] = nk + nj2 + i1m1 + 1;
			ipt[4] = nk2 + nj + i1m1;
			ipt[5] = nk2 + nj + i1m1 + 1;
			ipt[6] = nk2 + nj2 + i1m1;
			ipt[7] = nk2 + nj2 + i1m1 + 1;

			if (ido1d) {
				for (int iii = 0; iii < 4; iii++) {
					du[k1 - 1] += dlen * ds[iii];
					du[k1] += dlen * ds[iii + 4];
				}
			}
			else {
				for (int iii = 0; iii < 8; iii++)
					du[ipt[iii]] += ds[iii] * dlen;
			}

			if (iray) {
				double xray1 = xx / degrad;
				double yray1 = glatinv(hpi - yy) / degrad;
				double xray2 = xs / degrad;
				double yray2 = glatinv(hpi - ys) / degrad;
				FILE* fp_pth = fopen("", "");
				assert(0);
				fprintf(fp_pth, "%12.7f %12.7f %12.7f", xray1, yray1, zz);
				fprintf(fp_pth, "%12.7f %12.7f %12.7f", xray2, yray2, zs);
				fclose(fp_pth);
			}

			// ----rdevav is a measure of the deviation of ray from the source - station plane
			if (iraystat) {
				double rdevav = rdevs / nseg;
				rdevs = 0.;
				fprintf(fp_ray, " %10.4f %10.4f %10.4f %10.4f %10.4f %s %10.4f", ex, ey, ez, rsx, rsy, stt[kn], rdevav);
			}

			// ****** BEGIN STUFF FROM SPHYPIT
			if (nomat == 0) {
				//c------------------------------------------------------------------------------ -
				//c----- Comments on the meaning of variables
				//c  the following variables are accumulated for the observation
				//c     nbk = records the number of hits for a ray
				//c  the following are accumulated for each event
				//c     ind(i, j) = the block number of the block i hit by obs j
				//c     vm(i, j) = derivative of block i for obs j(weighted)
				//c     dat2(i) = temporary data vector
				//c     kbl = increment for unknowns, total number of blocks hit
				//c     mndex(i) = records position of new unknown
				//c     inbk(j) = number of blocks hit by obs j
				//c     mndm(i) = inverse of mndex, tells unknown in a position
				//c  the following are accumulated for the entire problem
				//c     nhit(i) = all hits in block i
				//c     indx(i) = position of new unknown
				//c     mbl = total unknowns
				//c----------------------------------------------------------------------------------
				if (i == 0) {
					kbl = 0;
					for (i1 = 0; i1 < nxyzc2; i1++) {
						mhit[i1] = 0;
					}
				}

				//c----Note on weighting : we weight each equation by sqrt(wt) because the
				//c      least squares solution will solve GT(Cdd) - 1G x = Gt(Cdd) - 1d.The elemenets
				//c       of Cdd are the variances(inverse is stored in pwt).We multiply both G and d by
				//c       sqrt(pwt) so that subsequent multipication recovers(Cdd) - 1.
				//
				//c---- - data vector : Note that in cases where we fit the deviation about
				//c   the mean we do not apply weights until after the residuals have been
				//c   demeaned.

				if ((isshot == 1 && idmean == 1) || istel == 1) {
					dat[i] = dt;
				}
				else {
					dat[i] = dt * wet;
				}

				//----- form hypocenter matrix
				am[i][0] = wet;
				am[i][1] = dtdx * wet;
				am[i][2] = dtdy * wet;
				am[i][3] = dtdz * wet;

				// ----- form velocity matrix
				int nbk = 0;
				for (int jb = 0; jb < maxvarc; jb++) {
					if (du[jb] != 0) {
						j = jb;

						// !if for S wave, add j, lzw
						if (iph == 1) {
							j += maxvarc;
						}

						// !===for Vp and Vs derivative, lzw
						vm[nbk][i] = du[jb] * wet;
						ind[nbk][i] = j;
						nbk++;
						if (nbk > maxnbk) {
							printf(" Error: nbk too large, stopping\n");
							fprintf(fp_log, " Error: nbk too large, stopping\n");
							assert(0);
						}

						if (mhit[j] == 0) {
							mndm[kbl] = j;
							mndex[j] = kbl;
							mhit[j] = 1;

							kbl++;
							if (kbl > maxkbl) {
								printf(" Error: kbl=%d\n", kbl);
								fprintf(fp_log, " Error: kbl=%d\n", kbl);
								assert(0);
							}
						}
						if (nhit[j] == 0) {
							indx[j] = mbl;
							jndx[mbl] = j;
							mbl++;
							if (mbl > maxmbl) {
								printf(" Error: mbl too large, stopping\n");
								fprintf(fp_log, " Error: mbl too large, stopping\n");
								assert(0);
							}
						}
						nhit[j]++;

						//c-- - Vp / Vs section : If we solve for dr instead of dUs, then
						//c	  Tso - Tsc = sum(dt / dUs) DUs = sum(dt / dUs) D(rUp) =
						//c	  	    sum(dt / dUs)[UpDr + rDUp].
						//c
						//c	  Thus, we scale the above by Up, and add a DUp scaled by r
						//c	  lzw   Only for the S wave raypath, and add some new factor to Vp part

						// if is the S wave
						if (ivpvs == 1 && iph == 1) {
							int jbp = jb;
							int jbs = j;
							if (ido1d) {
								jbp = jb * nxyc;
								jbs = jbp + nxyzc;
							}

							// scale Dr part using P velocity
							vm[nbk][i] = vm[nbk][i] * sp[jbp] * vpvsscale;


							nbk++;
							if (nbk > maxnbk) {
								printf(" Error: nbk too large, stopping\n");
								fprintf(fp_err, " Error: nbk too large, stopping\n");
								fprintf(fp_log, " Error: nbk too large, stopping\n");
								assert(0);
							}

							// ----Scale DUp derivative by r, and note that this is a new hit variable(Up)
							// 	            vm(nbk, i) = du(jb)*wet*sp(j)
							vm[nbk][i] = du[jb] * wet * sp[jbs];
							ind[nbk][i] = jb;

							if (mhit[jb] == 0) {
								mndm[kbl] = jb;
								mndex[jb] = kbl;
								mhit[jb] = 1;
								kbl++;
								if (kbl > maxkbl) {
									printf(" Error: kbl=%d\n", kbl);
									fprintf(fp_err, " Error: kbl=%d\n", kbl);
									fprintf(fp_log, " Error: kbl=%d\n", kbl);
									assert(0);
								}
							}

							if (nhit[jb] == 0) {
								jndx[jb] = mbl;
								jndx[mbl] = jb;
								mbl++;
								if (mbl > maxmbl) {
									printf(" Error: mbl too large, stopping\n");
									fprintf(fp_err, " Error: mbl=%d\n", mbl);
									fprintf(fp_log, " Error: mbl too large, stopping\n");
									assert(0);
								}
							}
							nhit[jb]++;
						} //end if (ivpvs && iph)
					} //end if(fabs(du)<0.0001)
				} //end jb for
				inbk[i] = nbk;
			} // ****** END STUFF FROM SPHYPIT
			goto a200; // ----go do next phase
		}

		// FIND THE RAY[GRAD(T)]
		nk = nxy * ke;
		nj = nx * je;
		nk2 = nx * ny * (ke + 1);
		nj2 = nx * (je + 1);

		ipt[0] = nk + nj + ie;
		ipt[1] = nk + nj + ie + 1;
		ipt[2] = nk + nj2 + ie;
		ipt[3] = nk + nj2 + ie + 1;
		ipt[4] = nk2 + nj + ie;
		ipt[5] = nk2 + nj + ie + 1;
		ipt[6] = nk2 + nj2 + ie;
		ipt[7] = nk2 + nj2 + ie + 1;
		//c  Discussion : The wavefront is a zero time change surface, which means
		//c	that the ray, perpendicular to the wavefront, is a maximum time
		//c	change direction.We therefore calculate the time gradient in the
		//c	cell to determine the local ray direction.
		//c
		//c  Note that because we are tracing the ray to the source, we are going
		//c	backwards in time, and hence need to look for the maximum negative
		//c	gradient.
		//c
		//c(1 / rsinQ)dt / df = [(pt2 - pt1) / (r1sinq1)+(pt4 - pt3) / (r1sinq2)
		//c + (pt6 - pt5) / (r2sinq1)+(pt8 - pt7) / (r2sinq2)] / 4df
		//c(1 / r)dt / dq = [(pt3 - pt1 + pt4 - pt2) / (r1)
		//c + (pt7 - pt5 + pt8 - pt6) / (r2)] / 4dq
		//c dt / dr = [(pt5 - pt1 + pt6 - pt2
		//	c + pt7 - pt3 + pt8 - pt4)] / 4h

		r1 = rearth - zk;
		r2 = r1 - h;
		q1 = yj;
		q2 = q1 + dq;
		f1 = xi;
		f2 = f1 + dq;

		// dt/df=(1/rsinQ)dt/df
		gradt[0] = ((t[j][ipt[1]] - t[j][ipt[0]]) / (r1 * sin(q1))
			+ (t[j][ipt[3]] - t[j][ipt[2]]) / (r1 * sin(q2))
			+ (t[j][ipt[5]] - t[j][ipt[4]]) / (r2 * sin(q1))
			+ (t[j][ipt[7]] - t[j][ipt[6]]) / (r2 * sin(q2))) / (4 * df);

		// dt / dq' = (1/r)dt/dq
		gradt[1] = ((t[j][ipt[2]] + t[j][ipt[3]]
			- t[j][ipt[0]] - t[j][ipt[1]]) / r1
			+ (t[j][ipt[6]] + t[j][ipt[7]]
				- t[j][ipt[4]] - t[j][ipt[5]]) / r2) / (4 * dq);

		// dt / dr' = -dt/dr
		gradt[2] = (t[j][ipt[4]] + t[j][ipt[5]]
			+ t[j][ipt[6]] + t[j][ipt[7]]
			- t[j][ipt[0]] - t[j][ipt[1]]
			- t[j][ipt[2]] - t[j][ipt[3]]) / (4 * h);

		//c-- - reverse sign on gradient since we are tracing the ray BACKWARDS(i.e.to the source, j).gradt3 is
		//c   a double negative(-dtdz = -(-dtdr)) so we leave it as is, but recall that gradt(3) is now the r component
		//c   of the ray.
		gradt[0] = gradt[0] * -1;
		gradt[1] = gradt[1] * -1;

		double sinf = sin(xx);
		double cosf = cos(xx);
		sinq = sin(yy);
		cosq = cos(yy);
		ro = rearth - zz;
		xo = ro * sinq * cosf;
		yo = ro * sinq * sinf;
		zo = ro * cosq;

		// if ray in seismometer cube, use straight ray from source
		if (ie >= (is - 2) && ie < (is + 2) &&
			je >= (js - 2) && je < (js + 2) &&
			ke >= (ks - 2) && ke < (ks + 2)) {
			r2 = rearth - zs;
			double x2 = r2 * sin(ys) * cos(xs);
			double y2 = r2 * sin(ys) * sin(xs);
			double z2 = r2 * cos(ys);
			double dx = x2 - xo;
			double dy = y2 - yo;
			double dz = z2 - zo;
			dl = sqrt(dx * dx + dy * dy + dz * dz);
			sx = dx / dl;
			sy = dy / dl;
			sz = dz / dl;
		}
		else {

			//c  NB : Because the spherical system is not orthoganal, the use of
			//c	Pathagoras here is not really legitimate.However, since we don't
			//c	use gradtm for anything after this step, all it does is scale the S
			//c	vector in spherical coordinates.We take the extra step of normalizing
			//c	S in the cartesian system by computing smag below.
			//c	gradtm is thus not really necessary, and is kept here only for reasons
			//c	of heritage.

			gradtm = sqrt(gradt[0] * gradt[0] + gradt[1] * gradt[1] + gradt[2] * gradt[2]);
			sf = gradt[0] / gradtm;
			sq = gradt[1] / gradtm;
			sr = gradt[2] / gradtm;
			// Rotate to the cartesian frame
			sx = -sinf * sf + cosq * cosf * sq + sinq * cosf * sr;
			sy = cosf * sf + cosq * sinf * sq + sinq * sinf * sr;
			sz = -sinq * sq + cosq * sr;
			smag = sqrt(sx * sx + sy * sy + sz * sz);
			sx /= smag;
			sy /= smag;
			sz /= smag;
		}

		//----save ray direction from source for use in focal mechanisms
		if (nseg == 0) {
			double sfq = sqrt(sf * sf + sq * sq);
			double az = atan2(sf, -sq) / degrad;
			double ai = atan2(sfq, -sr) / degrad;
			ais[i] = ai;
			azs[i] = az;
		}

		//----- intersection points for R constant surfaces
		dd[3];
		dd[2] = -1e20;
		if (gradt[2] != 0) {
			if (gradt[2] > 0)
				ro = rearth - zk;
			else
				ro = rearth - zk - h;
			b = xo * sx + yo * sy + zo * sz;
			c = xo * xo + yo * yo + zo * zo - ro * ro;
			quad = b * b - c;
			if (quad >= 0) {
				dd[2] = pmin(-b + sqrt(quad), -b - sqrt(quad));
			}
		}

		//----- intersection points for Q constant surfaces
		dd[1] = -1e20;
		if (gradt[1] != 0) {
			if (gradt[1] < 0) {
				sinq = sin(yj);
				cosq = cos(yj);
			}
			else {
				sinq = sin(yj + dq);
				cosq = cos(yj + dq);
			}
			a = 0;
			b = 0;
			c = 0;
			if (fabs(cosq) > fabs(sinq)) {
				tanq = sinq / cosq;
				tansq = tanq * tanq;
				a = sx * sx + sy * sy - sz * sz * tansq;
				b = xo * sx + yo * sy - zo * sz * tansq;
				c = xo * xo + yo * yo - zo * zo * tansq;
			}
			else {
				ctanq = cosq / sinq;
				ctansq = ctanq * ctanq;
				a = (sx * sx + sy * sy) * ctansq - sz * sz;
				b = (xo * sx + yo * sy) * ctansq - zo * sz;
				c = (xo * xo + yo * yo) * ctansq - zo * zo;
			}

			quad = b * b - a * c;
			if (quad >= 0) {
				dd[1] = pmin((-b + sqrt(quad)) / a, (-b - sqrt(quad)) / a);
			}
		}

		//------intersection points for F constant surfaces
		dd[0] = -1e20;
		if (gradt[0] != 0) {
			if (gradt[0] < 0) {
				sinf = sin(xi);
				cosf = cos(xi);
			}
			else {
				sinf = sin(xi + df);
				cosf = cos(xi + df);
			}
			double d1;
			if (fabs(cosf) > fabs(sinf)) {
				double tanf = sinf / cosf;
				d1 = (xo * tanf - yo) / (sy - sx * tanf);
			}
			else {
				double ctanf = cosf / sinf;
				d1 = (yo * ctanf - xo) / (sx - sy * ctanf);
			}
			if (d1 > 0) {
				dd[0] = d1;
			}
		}

		// ---find the possible step scalars
		if (gradt[0] < 0) {
			if (dd[0] < h / 1000) {
				//
				// ---The following tests prevent the ray from doubling back on itself.  The
				// 	  first is simply to keep the ray from going back to where it's just been.
				// 	  The second prevents a loop over 4 cells.   In this particular blockif,
				// 	  iseg will be -1, so a previous iseg of 1 means that the ray is simply
				// 	  trying to return to the cell it just left.  An example of a 4 cell loop
				// 	  would be isegs of 20, 1, -20, where an additional -1 would bring it
				// 	  back to the original cell, and hence set up the possibility of a closed
				// 	  loop.
				//
				if ((iseg == 1) || ((iseg + isegm1 + isegm2 == 1) && (nseg >= 3))) {
					gradt[0] = 0.;
					dd[0] = -1e20;
				}
			}
		}
		else if (gradt[0] > 0) {
			if (dd[0] < h / 1000) {
				if ((iseg == -1) || ((iseg + isegm1 + isegm2 == -1) && (nseg >= 3))) {
					gradt[0] = 0;
					dd[0] = -1e20;
				}
			}
		}
		else {
			dd[0] = -1e20;
		}

		// ################## 1
		if (gradt[1] < 0) {
			if (dd[1] < h / 1000) {
				if ((iseg == 20) || ((iseg + isegm1 + isegm2 == 20) && (nseg >= 3))) {
					gradt[1] = 0.;
					dd[1] = -1e20;
				}
			}
		}
		else if (gradt[1] > 0) {
			if (dd[1] < h / 1000) {
				if ((iseg == -20) || ((iseg + isegm1 + isegm2 == -20) && (nseg >= 3))) {
					gradt[1] = 0;
					dd[1] = -1e20;
				}
			}
		}
		else {
			dd[1] = -1e20;
		}

		// ################## 2
		if (gradt[2] > 0) {
			if (dd[2] < h / 1000) {
				if ((iseg == 300) || ((iseg + isegm1 + isegm2 == 300) && (nseg >= 3))) {
					gradt[2] = 0.;
					dd[2] = -1e20;
				}
			}
		}
		else if (gradt[2] < 0) {
			if (dd[2] < h / 1000) {
				if ((iseg == -300) || ((iseg + isegm1 + isegm2 == -300) && (nseg >= 3))) {
					gradt[2] = 0;
					dd[2] = -1e20;
				}
			}
		}
		else {
			dd[2] = -1e20;
		}

		// determine which is desired
		dfind(dd, &d, &md, tolmin, tolmax);

		dlen = d;
		xold = xx;
		yold = yy;
		zold = zz;

		xn = xo + d * sx;
		yn = yo + d * sy;
		zn = zo + d * sz;

		xx = atan2(yn, xn);
		yy = atan2(sqrt(xn * xn + yn * yn), zn);
		zz = rearth - sqrt(xn * xn + yn * yn + zn * zn);
		// printf("2034 xn=%lf\tyn=%lf\tzn=%lf\txx=%lf\tyy=%lf\tzz=%E\n", xn, yn, zn, xx, yy, zz);

		nseg++;
		isegm2 = isegm1;
		isegm1 = iseg;

		// ----x side intersected
		if (md == 0) {
			je = (int)(((yy - y[0]) / dq));
			ke = (int)(((zz - z0) / h));
			if (gradt[0] < 0) {
				if (ie < 0) {
					goto a150;
				}
				ie--;
				iseg = -1;
			}
			else {
				if (ie >= nx - 1) {
					goto a150;
				}
				ie++;
				iseg = 1;
			}
			// ----y side intersected
		}
		else if (md == 1) {
			if (gradt[1] < 0) {
				if (je < 0)
					goto a150;
				je--;
				iseg = -20;
			}
			else {
				if (je >= ny - 1) {
					goto a150;
				}
				je++;
				iseg = 20;
			}
			ie = (int)(((xx - x0) / df));
			ke = (int)(((zz - z0) / h));
			// ----z side intersected
		}
		else {
			if (gradt[2] > 0) {
				if (ke < 0) {
					goto a150;
				}
				ke--;
				iseg = -300;
			}
			else {
				if (ke >= nz - 1) {
					goto a150;
				}
				ke++;
				iseg = 300;
			}
			ie = (int)(((xx - x0) / df));
			je = (int)(((yy - y[0]) / dq));
		}
		if ((xx < x0 || ie >= nx) || (yy < y[0] || je >= ny) || (zz < z0 || ke >= nz)) {
			goto a150;
		}

		xm = (xx + xold) * 0.5;
		ym = (yy + yold) * 0.5;
		zm = (zz + zold) * 0.5;

		// printf("2103 xm=%lf\tym=%lf\tzm=%lf\n", xm, ym, zm);

		int i1 = -1, j1 = -1, k1 = -1;

		//--- find xm in coarse grid
		for (int iii = 1; iii < nxc; iii++) {
			if (xm < gx[iii]) {
				i1 = iii;
				break;
			}
		}
		if (i1 == -1) {
			printf(" Error xm = %lf is out of bounds!\n", xm);
			fprintf(fp_err, " Error xm = %lf is out of bounds!\n", xm);
			fprintf(fp_log, " Error xm = %lf is out of bounds!\n", xm);
			assert(0);
		}
		for (int iii = 1; iii < nyc; iii++) {
			if (ym < gy[iii]) {
				j1 = iii;
				break;
			}
		}
		if (j1 == -1) {
			printf(" Error ym = %lf is out of bounds!\n", ym);
			fprintf(fp_err, " Error ym = %lf is out of bounds!\n", ym);
			fprintf(fp_log, " Error ym = %lf is out of bounds!\n", ym);
			assert(0);
		}
		for (int iii = 1; iii < nzc; iii++) {
			if (zm < gz[iii]) {
				k1 = iii;
				break;
			}
		}
		if (k1 == -1) {
			printf(" Error zm = %lf is out of bounds!\n", zm);
			fprintf(fp_err, " Error zm = %lf is out of bounds!\n", zm);
			fprintf(fp_log, " Error zm = %lf is out of bounds!\n", zm);
			assert(0);
		}

		double hx = gx[i1] - gx[i1 - 1];
		double hy = gy[j1] - gy[j1 - 1];
		double hz = gz[k1] - gz[k1 - 1];
		float xci = gx[i1];
		float ycj = gy[j1];
		float zck = gz[k1];

		// --- derivatives
		fx = (xm - xci) / hx;
		fy = (ym - ycj) / hy;
		fz = (zm - zck) / hz;

		ds[0] = (1. - fz) * (1. - fy) * (1. - fx);
		ds[1] = (1. - fz) * (1. - fy) * fx;
		ds[2] = (1. - fz) * fy * (1. - fx);
		ds[3] = (1. - fz) * fy * fx;
		ds[4] = fz * (1. - fy) * (1. - fx);
		ds[5] = fz * (1. - fy) * fx;
		ds[6] = fz * fy * (1. - fx);
		ds[7] = fz * fy * fx;

		nk = nxyc * (k1 - 1);
		nj = nxc * (j1 - 1);
		nk2 = nxyc * k1;
		nj2 = nxc * j1;
		int i1m1 = i1;

		int jpt[8];
		jpt[0] = nk + nj + i1m1;
		jpt[1] = nk + nj + i1m1 + 1;
		jpt[2] = nk + nj2 + i1m1;
		jpt[3] = nk + nj2 + i1m1 + 1;
		jpt[4] = nk2 + nj + i1m1;
		jpt[5] = nk2 + nj + i1m1 + 1;
		jpt[6] = nk2 + nj2 + i1m1;
		jpt[7] = nk2 + nj2 + i1m1 + 1;
		if (ido1d == 0) {
			for (int iii = 0; iii < 8; iii++) {
				du[jpt[iii]] += ds[iii] * dlen;
			}
		}
		else {
			for (int iii = 0; iii < 4; iii++) {
				du[k1 - 1] += dlen * ds[iii];
				du[k1] += dlen * ds[iii + 4];
			}
		}

		if (iray == 1) {
			FILE* fp_pth = NULL;
			if (nseg == 1) {
				if (isshot == 0) {
					if (istel == 0) {
						sprintf(rayfile, "%s_event%d_%c.ray\n", sta[i], nev, phs[i]);
					}
					else {
						sprintf(rayfile, "%s_tele%d_%c.ray\n", sta[i], ntel, phs[i]);
					}
				}
				else {
					sprintf(rayfile, "%s_shot%d_%c.ray\n", sta[i], nshot, phs[i]);
				}
				d_blank(rayfile, &len);
				printf("File: %s\n", rayfile);

				fp_pth = fopen(rayfile, "r");
				if (!fp_pth) {
					printf("error on opening file (%s)\n", rayfile);
					assert(0);
				}
			}
			double xray1 = xold / degrad;
			double yray1 = glatinv(hpi - yold) / degrad;
			double xray2 = xx / degrad;
			double yray2 = glatinv(hpi - yy) / degrad;
			if (fp_pth) {
				fprintf(fp_pth, "%12.5f %12.5f %12.5f\n", xray1, yray1, zold);
				fprintf(fp_pth, "%12.5f %12.5f %12.5f\n", xray2, yray2, zz);
			}
		}

		// ---temporarly output for piercing point figure
		// 	  write (31,fmt='(3f10.3)') xold,yold,-zold

		// ----rdev is a measure of the deviation of ray from the source-station plane
		double rx = xm - ex;
		double ry = ym - ey;
		double rmagsq = rx * rx + ry * ry;
		double rds = rsx * rx + rsy * ry;
		double rdev = sqrt(rmagsq - rds * rds);
		rdevs += rdev;
		// -----------------------------------------
		// ---go back to get next ray segment
		goto a100;
	a150:
		printf(" Warning: Ray has left the medium.  Skipping remainder\n");
		fprintf(fp_err, " Warning: Ray has left the medium.  Skipping remainder\n");
		isgood[i] = 0;
	a200: {
		}
	} // ***** End of Loop over Observations *****

	// --- demean the data vector if required
	avresev /= facsev;
	if ((isshot == 1 && idmean == 1) || istel == 1) {
		if (isshot == 1) {
			printf(" Average residual for shot %d is %lf", nshot, avresev);
		}
		if (istel == 1) {
			printf(" Average residual for tele %d is %lf", ntel, avresev);
		}
		dpot += avresev;
		for (k = 0; k < nsta; k++) {
			if (isgood[k] == 1) {
				dat[k] -= avresev;

				// --note multipyling by pwt will make sum of residuals be zero.
				//  However, in order to be consistent with the local / shot weighting, we
				//  multiply the data vectory by sqrt(pwt(k)) anticipating that another sqrt(wt)
				//  will come from a weighted GT.
				dat[k] = dat[k] * sqrt(pwt[k]);
				resmin[k] = resmin[k] - avresev;
				obstime[k] = obstime[k] - avresev;
			}
		}
	}

	// *** BEGIN STUFF FROM SPHYPIT
	if (nomat == 0) {
		for (int j1 = 0; j1 < nsta; j1++) {
			for (int j2 = 0; j2 < kbl; j2++) {
				vmp[j2][j1] = 0;
			}
		}

		// Note that for shot data, if we demean the residuals we should also demean the
		// partial derivatives, just like with teleseisms.
		if ((isshot == 1 && idmean == 1) || istel == 1) {
			for (k = 0; k < kbl; k++) {
				vsum[k] = 0;
			}
			wsum = 0;

			// --Note on weighting:  we ultimately want vmp to be
			//
			//        vmp = sqrt(wt) * [ g - sum(wt*g)/sum(wt)]
			//
			//        Again, the sqrt(wt) at the beginning is because of postmultiplication in
			//        Least Squares to construct (Cdd)-1 properly.  The data vector is multiplied
			//        above by sqrt(wt) as well.
			//
			//        At this point, vm is sqrt(wt)*g.  vsum will wind up
			//        being the average, so we sum sqrt(wt)*vm to get wt*g,
			// 	vmp is then sqrt(wt)*sqrt(wt)*g = wt*g
			//        At the end of this loop we have
			//        vsum(nvar) = sum over nobs (wt(nobs)*g(nvar, nobs))
			//        vmp(nvar, nobs) =  wt(nobs)*g(nvar,nobs)
			//        and wsum is just sum over nobs (wt(nobs)).
			for (k = 0; k < nsta; k++) {
				if (isgood[k]) {
					float wt = pwt[k];
					double wet = sqrt(wt);
					int lim = inbk[k];
					if (lim != 0) {
						for (int j1 = 0; j1 < lim; j1++) {
							int jnd = mndex[ind[j1][k]];
							float valv = vm[j1][k];
							vmp[jnd][k] = valv;
							// vsum = sum over j of Wij*gijl
							vsum[jnd] += valv * wet;
						}
						wsum += wt;
					}
				}
			}

			// At this point, vmp is wt*g, vsum is sum(wt*g), and wsum = sum(wt)
			for (int jnd = 0; jnd < kbl; jnd++) {
				//--- valv is now sum(wt*g)/sum(wt)
				float valv = vsum[jnd] / wsum;
				for (int i1 = 0; i1 < nsta; i1++) {
					if (isgood[i1]) {
						vmp[jnd][i1] -= sqrt(pwt[i1]) * valv;
					}
				}
			}
		}
		else {
			for (k = 0; k < nsta; k++) {
				if (isgood[k]) {
					int lim = inbk[k];
					if (lim != 0) {
						for (int j1 = 0; j1 < lim; j1++) {
							int jnd = mndex[ind[j1][k]];
							vmp[jnd][k] += vm[j1][k];
						}
					}
				}
			}
		}

		// ----station corrections terms
		//   Some explanation:
		//     We consider that for station j there is a True station correction
		//     Cj and a current esimate of it cj. We define a station correction
		//     as a number added to the ray travel time to obtain the "real" time.
		//     This is why the time correction is added above (dt = dt + tc).
		//     Thus, a residual (Observed - calculated) will involve a difference
		//     between the real (or "observered") and estimated correction as
		// 	dCj = Cj - cj
		//     For teleseisms, we demean the residuals for a given event i,
		//     and for the station correction term this looks like:
		//
		// 	Rij - 1/N(Ej Rij) = dCj - 1/N(Ej dCj)
		//
		//     In other words, we demean the correction term as well. In the
		//     normal equations, these terms look like
		//
		// 	0  1-1/N 0 0 -1/N  -1/N 0 0
		//
		//     The column corresponding to the observation is 1-1/N and the remaining
		//     observations of this event are -1/N.
		//
		//     When weighting, we replace N by wsum and the "1s" by the square
		//     root of individual weights (as above with vmp).
		//
		// 	NB: this part has not been properly tested.

		if (istel == 1 && istacor == 1) {
			ncwrt = 0;
			for (int k1 = 0; k1 < nsta; k1++) {
				for (int j1 = 0; j1 < nstr2; j1++) {
					dtc[j1] = 0;
				}
				if (isgood[k1]) {
					int ja = 0;
					for (int j1 = 0; j1 < nsta; j1++) {
						if (isgood[j1]) {
							jsave[ja] = j1;
							ja++;
							double sqrpwt = sqrt(pwt[j1]);
							dtc[ja] = -sqrpwt / wsum;
							if (k1 == j1) {
								dtc[ja] = sqrpwt + dtc[ja];
							}
						}
					}
					fprintf(fp_stc, "%d\n", ja);
					for (int i1 = 0; i1 < ja; i1++) {
						fprintf(fp_stc, "%d %f\n", istn[jsave[i1]], dtc[i1]);
					}
					ncwrt++;
				}
			}
		}
		// Note: need equivalent to above for local events eventually.

		// ----build equations for path terms
		if (isshot == 0) {
			if (istel == 0) {
				printf(" Making path terms for eq = %d\n", nev);
			}
			else {
				printf("Making path terms for tele = %d\n", ntel);
			}
		}
		else {
			printf("Making path terms for shot = %d\n", nshot);
		}

		int nwrt = 0;
		for (int j1 = 0; j1 < nsta; j1++) {
			if (isgood[j1]) {
				int ja = 0;
				for (int i1 = 0; i1 < kbl; i1++) {
					float vm2 = vmp[i1][j1];
					if (vm2 != 0) {
						jsave[ja++] = i1;
					}
				}
				if (ja > 0) {
					fwrite(&ja, sizeof(ja), 1, fp_dts);
					for (int i1 = 0; i1 < ja; i1++) {
						indx[mndm[jsave[i1]]]++;
						fwrite(&indx[mndm[jsave[i1]]], sizeof(indx[mndm[jsave[i1]]]), 1, fp_dts);
						indx[mndm[jsave[i1]]]--;
						fwrite(&vmp[jsave[i1]][j1], sizeof(vmp[jsave[i1]][j1]), 1, fp_dts);
					}

					// write residual file
					fwrite(&dat[j1], sizeof(dat[j1]), 1, fp_res);

					if (isshot == 0) {
						if (istel == 0) {
							for (int i1 = 0; i1 < 4; i1++) {
								fprintf(fp_dth, "%f\n", am[j1][i1]);
							}
						}
						// ---in the case of a shot, we choose to make the origin time a variable
					}
					else {
						fprintf(fp_dth, "%f\n", am[j1][0]);
					}
					nwrt++;
				}
				else {
					printf(" Warning: ray has no derivative!\n");
					printf(" Potential problem for syncing with corrections!\n");
					fprintf(fp_err, " Warning: ray has no derivative!\n");
					fprintf(fp_err, " Potential problem for syncing with corrections!\n");
				}
			}
		}
		if (isshot == 0) {
			if (istel == 0) {
				fprintf(fp_bok, "%d %d\n", nev, nwrt);
			}
			else {
				fprintf(fp_bok, "%d %d\n", m2, nwrt);
				if (istacor && ncwrt != nwrt) {
					printf(" Error: ncwrt, nwrt not equal = %d %d\n", ncwrt, nwrt);
					fprintf(fp_err, " Error: ncwrt, nwrt not equal = %d %d\n", ncwrt, nwrt);
				}
			}
		}
		else {
			// ---in this case we solve for shot origin times
			// fprintf(fp_bok, "%d %d\n", m1, nwrt);
			// ---in this case we do not solve for shot origin times
			fprintf(fp_bok, "%d %d\n", m2, nwrt);
		}
	} // END STUFF FROM SPHYPIT

	// ----output a data file if requested
	if (idatout == 1) {
		int iyr = 0, jday = 0, ihr = 0, imn = 0;
		double sec = 0;
		etoh(dpot, &iyr, &jday, &ihr, &imn, &dsec);
		sec = dsec;
		if (istel == 0) {
			xlon = ex / degrad;
			xlat = glatinv(hpi - ey) / degrad;
		}
		fprintf(fp_hed, "%4d %3d %2d %2d %8.4lf %9.5f %10.5lf %8.4lf %12s %13.3lf\n", iyr, jday, ihr, imn, sec, xlat, xlon, ez, evid, avresev);
		fprintf(fp_dot, "%4d %3d %2d %2d %8.4lf %9.5f %10.5lf %8.4lf %12s %13.3lf\n", iyr, jday, ihr, imn, sec, xlat, xlon, ez, evid, avresev);
		for (j = 0; j < nsta; j++) {
			if (isgood[j]) {
				tarr = obstime[j] + dpot;
				etoh(tarr, &iyr, &jday, &ihr, &imn, &dsec);
				// ---call attention to residuals over resflag seconds in magnitude
				if (fabs(resmin[j]) > resflag) {
					mark = '*';
				}
				else {
					mark = ' ';
				}
				if (istel == 1) {
					fprintf(fp_dot, "%s  %4d %3d %2d %2d %7.3lf %c %15.3f %10.3lf %c %lf %lf %lf\n", sta[j], iyr, jday, ihr, imn, dsec, phs[j], rwts[j], resmin[j], mark, tdelts[j], az1s[j], az2s[j]);
				}
				else {
					fprintf(fp_dot, "%s  %4d %3d %2d %2d %7.3lf %c %15.3f %10.3lf %c %lf %lf\n", sta[j], iyr, jday, ihr, imn, dsec, phs[j], rwts[j], resmin[j], mark, azs[j], ais[j]);
				}
			}
		}
		fprintf(fp_dot, "\n");
	}

	goto a3; // 2474 go back and start on next earthquake
a60:
	{
		FILE *fp_sht = NULL, *fp_tel = NULL;
		// ---go back and do shots if necessary
		if (isshot == 0) {
			if (istel == 0) {
				istel = idotel;
				if (istel == 1) {
					fp_tel = fopen(telefil, "r");
					if (!fp_tel) {
						printf("error on opening file (%s)\n", telefil);
						assert(0);
					}
					fp_din = fp_tel;

					// ----- read in time tables - hardwire for test
					FILE* fp_ptm = fopen(pbasfil, "r");
					if (!fp_ptm) {
						printf("error on opening file (%s)\n", pbasfil);
						assert(0);
					}
					FILE* fp_stm = fopen(sbasfil, "r");
					if (!fp_stm) {
						printf("error on opening file (%s)\n", sbasfil);
						assert(0);
					}
					// readtbl3(fp_ptm, fp_stm);

					// ----- read in ellipticity correction table
					FILE* fp_elp = fopen(elipfil, "r");
					if (!fp_elp) {
						printf("error on opening file (%s)\n", elipfil);
						assert(0);
					}
					// redtab(fp_elp);

					// -----translate local grid coordinates to lat-lon for points at the base of the model.
					for (j = 0; j < ny; j++) {
						xlat = y[0] + dq * j;
						for (i = 0; i < nx; i++) {
							xlon = x0 + df * i;
							slons[i][j] = xlon;
							slats[i][j] = xlat;
						}
					}
					goto a3;
				}
				else {
					isshot = idoshot;
					if (isshot) {
						fp_sht = fopen(shotfil, "r");
						if (!fp_sht) {
							printf("error on opening file (%s)\n", shotfil);
							assert(0);
						}
						fp_din = fp_sht;
						istel = 0;
						goto a3;
					}
				}
			}
			else {
				isshot = idoshot;
				if (isshot) {
					fp_sht = fopen(shotfil, "r");
					if (!fp_sht) {
						printf("error on opening file (%s)\n", shotfil);
						assert(0);
					}
					fp_din = fp_sht;
					istel = 0;
					goto a3;
				}
			}
		}
		if (fp_eqs) {
			fclose(fp_eqs);
		}
		if (fp_sht) {
			fclose(fp_sht);
		}
		if (fp_tel) {
			fclose(fp_tel);
		}

		// ---- global statistics
		double avwt = facs / knobs;
		double deb = (rsq - avres * avres / facs) / (knobs * avwt);
		double std = sqrt(deb);
		rms = sqrt(rms / wtmod);
		avres = avres / facs;

		char form[34 * 6 + 31 * 4 + 1] = " Overall data variance :  %10.4f\n";
		strcat(form, "    standard deviation :  %10.4f\n");
		strcat(form, "                   RMS :  %10.4f\n");
		strcat(form, "      Average Residual :  %10.4f\n");
		strcat(form, "        Average Weight :  %10.4f\n");
		strcat(form, "            Chi-square :  %10.4f\n");
		strcat(form, "Number of Observations :  %6d\n");
		strcat(form, " Number of Earthquakes :  %6d\n");
		strcat(form, "  Number of Teleseisms :  %6d\n");
		strcat(form, "       Number of Shots :  %6d\n");
		printf(form, deb, std, rms, avres, avwt, sqrt(rsq) / knobs, knobs, nev, ntel, nshot);
		fprintf(fp_err, form, deb, std, rms, avres, avwt, sqrt(rsq) / knobs, knobs, nev, ntel, nshot);
		fprintf(fp_log, form, deb, std, rms, avres, avwt, sqrt(rsq) / knobs, knobs, nev, ntel, nshot);

		printf("\n");
		fprintf(fp_err, "\n");
		fprintf(fp_log, "\n");
		fprintf(fp_err, "  No. Stn    Phs Average St.Dev. Nobs\n");
		fprintf(fp_log, "  No. Stn    Phs Average St.Dev. Nobs\n");

		for (kn = 0; kn < nstr; kn++) {
			for (iph = 0; iph < 2; iph++) {
				facs = facstn[kn][iph];
				int jnobs = nobstn[kn][iph];
				rsq = rstnsq[kn][iph];
				avres = avrstn[kn][iph];
				if (facs > 0) {
					avwt = facs / knobs;
					deb = (rsq - avres * avres / facs) / (jnobs * avwt);
					std = sqrt(fabs(deb));
					avres = avres / facs;
				}
				else {
					avres = 0.;
					std = 0.;
				}
				fprintf(fp_err, "%d %s %d %lf %lf %d\n", kn, stt[kn], iph, avres, std, jnobs);
				fprintf(fp_log, "%d %s %d %lf %lf %d\n", kn, stt[kn], iph, avres, std, jnobs);
			}
		}

		// ----output indexing array
		if (nomat == 0) {
			// -----tack on a marker to signify end of file
			int ja = -1;
			fwrite(&ja, sizeof(ja), 1, fp_dts);
			fwrite(&mbl, sizeof(mbl), 1, fp_dts);
			for (int i=0;i < mbl ;i++){
				fwrite(&i, sizeof(i), 1, fp_dts);
				fwrite(&jndx[i], sizeof(jndx[i]), 1, fp_dts);
			}

			fprintf(fp_msc, "%d\n", mbl);
			for (i = 0; i < mbl; i++) {
				fprintf(fp_msc, "%d\n", jndx[i]);
			}
			if (istacor) {
				fprintf(fp_stc, "%d\n", ja);
			}
			printf("\n");
			printf(" Number of Wavespeed Variables: %d\n", mbl);
			printf(" Degrees of freedom: %d\n", knobs - mbl);
			fprintf(fp_err, "\n");
			fprintf(fp_err, " Number of Wavespeed Variables: %d\n", mbl);
			fprintf(fp_err, " Degrees of freedom: %d\n", knobs - mbl);

			// ----ouput nhit array
			fprintf(fp_hit, "%d", mbl);
			for (i = 0; i < mbl; i++) {
				fprintf(fp_hit, "%d %d\n", jndx[i], nhit[jndx[i]]);
			}

			// ----put a trailer on the hypo bookeeping file to signify the end of data
			// ---we use m2 because the hypo programs will presume this is a tele.
			fprintf(fp_bok, "%d %d\n", m2, m2);
		}

		fclose(fp_dts);
		fclose(fp_msc);
		fclose(fp_err);
	}
	return 0;
}
