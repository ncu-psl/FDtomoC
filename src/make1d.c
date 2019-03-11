/*
c---General code to generate a coarse mesh 1D model from an
c	ASCII file specification and a spec file.
c
c	Author:	Steve Roecker, RPI
c
c	Version 2004.0909
c
c	Input files:
c
c	 Name   Uname   U#  Comments
c       onedfil lunone  14  ASCII file describing 1D model
c
c	Output files:
c
c	 Name   Uname   U#  Comments
c       oldvfil luncor  12  Coarse grid wavespeed model file
c
c	oldvfil is specified in the Spec file.
c
c	The format of the  1D input file should be:
c
c	<Descriptive Header>
c	Z1	Vp1	Vs1	C/I
c	Z2	Vp2	Vs2	C/I
c	...
c
c	where
c	Vp	is the P wavespeed
c	Vs	is the S wavespeed (or Vp/Vs if vs1d = 0)
c	Z	is the starting (upper) depth
c	C/I	is the interpolation specifier.  Set
c		to C for constant speed to next depth, and to
c		I to interpolate to the next depth.  If
c		the last wavespeed in the list is set to I, then
c		speeds at deeper depths are extrapolated using the
c		last available gradient (the difference between the
c		last two speeds divided by the difference between the
c		last two depths).
c
c	Variables required by this program:
c
c	nxc	Number of coarse grid points in the x direction
c	nyc	Number of coarse grid points in the y direction
c	nzc	Number of coarse grid points in the z direction
c	h	Find grid spacing
c
c	Optional variables that can be usefully set by the user:
c
c	igridx	Array of spaces between fine and coarse grid in x direction
c	igridy	Array of spaces between fine and coarse grid in y direction
c	igridz	Array of spaces between fine and coarse grid in z direction
c	flat	= 1 to flatten depths and speeds, 0 to leave as is. NB:  Because the tomo
c			sofware uses a regular grid, we presume that the input grid description 
c			is already flattened.  This program starts by "unflattening" these
c			depths prior to looking up the appropriate wavespeed from the table.
c			The output is flattened wavespeeds.
c	vs1d	= 1 to interpret Vs as shear wave speed; = 0 to interpret as Vp/Vs
c	isph	= 1 to use a spherical coordinate system; =0 for cartesian.
c		    note that the only practical difference is that the definitions of dx and dy
c		    change so that dx = df and dy = df.
*/

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "common/environment_setting.h"
#include "common/parseprogs.h"
#include "common/parameter.h"
#include "common/gridspec.h"
#include "common/geographic_method.h"
#include "common/string_process.h"

#define MAX1D 1000
#define MAXSTRLEN 132

//--mustv is the number of variables that must be assigned in the parameter
//      file in order for this program to run correctly.  The names of the
//      variables are specified in the mvals string array below.
//--mustf is the number of files that must be attached; see the files array below.
#define MUSTV  4
#define MUSTF  2

//---number of 4 byte words in the header
#define nhbyte 58 * 4

double gx[nxcm], gy[nxcm], gz[nxcm];
double vp[MAX1D][2], z[MAX1D];

int iflat = 0, isph = 0;
double z0r;
double hpi = 1.570796;
double degrad = (double) 0.017453292;
double y00;

float vsave[nxyzcm2];
int vs1d = 1;

char spec_file[MAXSTRLEN + 1];
char aline[MAXSTRLEN + 1], varname[MAXSTRLEN + 1], parval[MAXSTRLEN + 1];
char *mvals[MUSTV] = { "nxc\0", "nyc\0", "nzc\0", "h\0" };
char *files[MUSTF] = { "oldvfil\0", "onedfil\0" };

char oldvfil[MAXSTRLEN + 1], onedfil[MAXSTRLEN + 1];

char VERSION[10] = "2017.1122\0";
char terp[MAX1D + 1];

//----header stuff
char head[5], type[5], syst[5];
char quant[5];
char flatten[5];
char hcomm[101];

//----fxs, fys, and fzs are not used in wavespeed models, so just set to zero
double fxs = 0.0, fys = 0.0, fzs = 0.0;
double clat, clon, cz;
double axo, ayo, azo, dx, dy, dz;
float az;
int nxh, nyh, nzh;
//--------------------------------------
char hdr[nhbyte + 1];

double rearth = 6371.0;
int lenhead = nhbyte * 4;

FILE *fp_log;
FILE *fp_spc;
FILE *fp_cor;
FILE *fp_one;

double flatvel(double, double);
double uflatz(double);
double flatz(double);
char * dtoa(char *, double, int);

int main() {
	char pval[MAXSTRLEN + 1];
	int len, ierr;
	printf("Enter parameter specification file: ");
	scanf("%s",spec_file);
	spec_file[MAXSTRLEN] = '\0';

	fp_spc = fopen(spec_file, "r");
	int i;

//---recover the variables needed to run this program
//
//       nxc, nyc, nzc      coarse dimensions of the fine mesh used in the trt tables
//       h                  fine grid spacing
	for (i = 0; i < MUSTV; i++) {
		get_vars(fp_spc, mvals[i], pval, &len, &ierr);
		if (ierr == 1) {
			goto a50;
		}
		if (i == 0) {
			sscanf(pval, "%d", &nxc);
		} else if (i == 1) {
			sscanf(pval, "%d", &nyc);
		} else if (i == 2) {
			sscanf(pval, "%d", &nzc);
		} else if (i == 3) {
			sscanf(pval, "%lf", &h);
		}
	}

//----dimension check
	if (nxc > nxcm) {
		printf("nxc is too large, maximum is: %d\n", nxcm);
		assert(nxc <= nxcm);
	} else if (nyc > nycm) {
		printf("nyc is too large, maximum is: %d\n", nycm);
		assert(nyc <= nycm);
	} else if (nzc > nzcm) {
		printf("nzc is too large, maximum is: %d\n", nzcm);
		assert(nzc <= nzcm);
	}

	int lenf1, lenf2;
	for (i = 0; i < MUSTF; i++) {
		get_vars(fp_spc, files[i], pval, &len, &ierr);
		if (ierr == 1) {
			goto a51;
		}
		if (i == 0) {
			sscanf(pval, "%s", oldvfil);
			lenf1 = len;
		} else if (i == 1) {
			sscanf(pval, "%s", onedfil);
			lenf2 = len;
		}
	}

//--Optionally read in some variables
//---Coordinate origin (used in header)
	get_vars(fp_spc, "x0 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &x0);
	get_vars(fp_spc, "y0 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &y[0]);
	get_vars(fp_spc, "z0 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &z0);
	get_vars(fp_spc, "clat ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &clat);
	get_vars(fp_spc, "clon ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &clon);
	get_vars(fp_spc, "cz ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &cz);
	get_vars(fp_spc, "azmod ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &az);
	get_vars(fp_spc, "df ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &df);
	get_vars(fp_spc, "dq ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &dq);

//----flatness, Vs, and sph  flags
	get_vars(fp_spc, "flat ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &iflat);

	get_vars(fp_spc, "vs1d ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &vs1d);

	get_vars(fp_spc, "sph ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &isph);

//-----Grid specs
	int ib = 0, ie = 0, lenv = 0, nvl = 0;
	rewind(fp_spc);
	a11: get_line(fp_spc, aline, &ierr);
	aline[MAXSTRLEN] = '\0';
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
	int k;
	for (k = 1; k < nxc; k++) {
		ib = ie;
		get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
		sscanf(parval, "%d", &igridx[k]);
	}
	a12: rewind(fp_spc);
	a13: get_line(fp_spc, aline, &ierr);
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
	a14: rewind(fp_spc);

	a15: get_line(fp_spc, aline, &ierr);
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

//-----end of optional parameters

	int nxyc = nxc * nyc;
	int nxyzc = nxyc * nzc;
	int nxyzc2 = nxyzc * 2;

	if (isph == 1) {
		y00 = y[0] * degrad;
//		y00 = hpi - glat(y00)
		y00 = hpi - glath(y00, z0, &z0r);

//	---If dq and df have not been specified, { make them so that the
//	   interval at the surface is equal to h
		if (dq == 0.0)
			dq = h / rearth;
		if (df == 0.0)
			df = fabs(h / (rearth * sin(y00)));
		dy = dq / degrad;
		dx = df / degrad;
		printf("dx=%lf dy=%lf df=%lf dq=%lf\n", dx, dy, df, dq);
	} else {
		dx = h;
		dy = h;
	}
	nx = 1;
	ny = 1;
	nz = 1;
	gx[0] = x0;
	gy[0] = y[0];
	gz[0] = z0;

	for (int i = 1; i < nxc; i++) {
		nx = nx + igridx[i - 1];
		gx[i] = gx[i - 1] + dx * igridx[i - 1];
	}

	for (int i = 1; i < nyc; i++) {
		ny = ny + igridy[i - 1];
		gy[i] = gy[i - 1] + dy * igridy[i - 1];
	}

	for (int i = 1; i < nzc; i++) {
		nz = nz + igridz[i - 1];
		gz[i] = gz[i - 1] + h * igridz[i - 1];
	}

	int nxy = nx * ny;
	int nxyz = nxy * nz;

//	----dimension check
	if (nx > nxm) {
		printf(" nx is too large, maximum is: %d", nxm);
		assert(!(nx > nxm));
	}
	if (ny > nym) {
		printf(" ny is too large, maximum is: %d", nym);
		assert(!(ny > nym));
	}
	if (nz > nzm) {
		printf(" nz is too large, maximum is: %d", nzm);
		assert(!(nz > nzm));
	}

	int lengrd = 4 * (nxc + nyc + nzc - 3);
	int lenrec = lenhead + lengrd + 4 * nxyzc2;

	fp_log = fopen("../../make1d.log", "w");
	if(!fp_log) {
		printf("(Error in make1d.c)create fp_log file error.\n");
		assert(fp_log);
	}

	fprintf(fp_log, "  \n");
	fprintf(fp_log,
			" *************************************************************** \n");

	fprintf(fp_log, "          Parameters Set For This Run of make1d.f\n");

	fprintf(fp_log, "  \n");
	fprintf(fp_log, " VERSION: %s\n", VERSION);
	fprintf(fp_log, "  \n");
	fprintf(fp_log, " Current parameter specification file: %-40s\n",
			spec_file);
	fprintf(fp_log, "  \n");
	char out_str[MAXSTRLEN];
	dtoa(out_str, clat, 18);
	fprintf(fp_log, "  Latitude origin  (clat):   %s     \n", out_str);
	dtoa(out_str, clon, 18);
	fprintf(fp_log, "  Longitude origin (clon):   %s     \n", out_str);
	dtoa(out_str, cz, 18);
	fprintf(fp_log, "  Depth of  origin (cz)  :   %s     \n", out_str);
	dtoa(out_str, az, 10);
	fprintf(fp_log, "  Clockwise rotation (az):   %s    \n", out_str);
	fprintf(fp_log, "  \n");
	dtoa(out_str, x0, 18);
	fprintf(fp_log, "  Cartesian X origin (x0):   %s     \n", out_str);
	dtoa(out_str, y[0], 18);
	fprintf(fp_log, "  Cartesian Y origin (y0):   %s     \n", out_str);
	dtoa(out_str, z0, 18);
	fprintf(fp_log, "  Cartesian Z origin (z0):   %s     \n", out_str);
	if (isph == 0) {
		fprintf(fp_log, "  Coordinate system is CARTESIAN \n");
		fprintf(fp_log, "  Fine grid spacing: %lf\n", h);
	} else {
		fprintf(fp_log, "  Coordinate system is SPHERICAL \n");
		fprintf(fp_log, "  Fine Longitude spacing (df):    %.16E\n", dx);
		fprintf(fp_log, "  Fine Latidtude spacing (dq):    %.16E\n", dy);
		dtoa(out_str, h, 18);
		fprintf(fp_log, "  Fine Radial spacing    (dz):   %s     \n", out_str);
	}
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of X coarse grid nodes: %12d\n", nxc);
	fprintf(fp_log, "  X coarse grid node spacing: \n");
	for (int iii = 0; iii < nxc - 1; iii++) {
		fprintf(fp_log, "% 4d", igridx[iii]);
		if (iii % 10 == 9) {
			fprintf(fp_log, "\n");
		}
	}
	fprintf(fp_log, "\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of Y coarse grid nodes: %12d\n", nyc);
	fprintf(fp_log, "  Y coarse grid node spacing: \n");
	for (int iii = 0; iii < nyc - 1; iii++) {
		fprintf(fp_log, "% 4d", igridy[iii]);
		if (iii % 10 == 9) {
			fprintf(fp_log, "\n");
		}
	}
	fprintf(fp_log, "\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of Z coarse grid nodes: %12d\n", nzc);
	fprintf(fp_log, "  Z coarse grid node spacing: \n");
	for (int iii = 0; iii < nzc - 1; iii++) {
		fprintf(fp_log, "% 4d", igridz[iii]);
		if (iii % 10 == 9) {
			fprintf(fp_log, "\n");
		}
	}
	fprintf(fp_log, "\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of X fine grid nodes: %12d\n", nx);
	fprintf(fp_log, "  Number of Y fine grid nodes: %12d\n", ny);
	fprintf(fp_log, "  Number of Z fine grid nodes: %12d\n", nz);
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Total Number of coarse grid nodes:%13d\n", nxyzc);
	fprintf(fp_log, "  Total Number of fine grid nodes:%13d\n", nxyz);
	fprintf(fp_log, " \n");
	if (iflat == 1) {
		fprintf(fp_log, "  Speeds and Depths are flattened (iflat = 1) \n");
	} else {
		fprintf(fp_log, "  Speeds and Depths are not flattened (iflat = 0) \n");
	}
	if (vs1d == 1) {
		fprintf(fp_log, "  Vs column treated as S wavespeed (vs1d = 1) \n");
	} else {
		fprintf(fp_log, "  Vs column treated as Vp/Vs (vs1d = 0) \n");
	}

	fprintf(fp_log, " \n");

	fprintf(fp_log, "  Length of header  (bytes): %12d\n", lenhead + lengrd);
	fprintf(fp_log, "  Length of outfile (bytes): %12d\n", lenrec);

	fprintf(fp_log, " \n");
	fprintf(fp_log, " Input file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " One-D ASCII model file: %-40s\n", onedfil);
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Output file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " New Wavespeed model file: %-40s\n", oldvfil);
	fprintf(fp_log, " \n");

//---generate the header
	strcpy(head, "HEAD");
	head[4] = '\0';
	strcpy(type, "CORS");
	type[4] = '\0';
	if (isph == 1) {
		strcpy(syst, "SPHR");
		dx = df / degrad;
		dy = dq / degrad;
	} else {
		strcpy(syst, "CART");
		dx = h;
		dy = h;
	}
	syst[4] = '\0';
	strcpy(quant, "BMOD");
	quant[4] = '\0';
	if (iflat == 1) {
		strcpy(flatten, "FLAT");
	} else {
		strcpy(flatten, "NOFL");
	}
	flatten[4] = '\0';
	sprintf(hcomm, "Output from make1d.c using %s", onedfil);
	hcomm[100] = '\0';
	axo = x0;
	ayo = y[0];
	azo = z0;
	dz = h;
	nxh = nxc;
	nyh = nyc;
	nzh = nzc;

//---unflatten the depths if required
	if (iflat == 1) {
		for (i = 0; i < nzc; i++) {
			gz[i] = uflatz(gz[i]);
		}
	}

	fp_one = fopen(onedfil, "r");
	if(!fp_one) {
        printf("Error on opening fp_one(%s)\n", onedfil);
        assert(0);
	}
//---skip over header
	get_line(fp_one, aline, &ierr);
	int nl = 0;
	a21: get_line(fp_one, aline, &ierr);
	if (ierr == 1)
		goto a2;
	if (ierr != 0)
		goto a21;
	ib = 0;
	//******************************************
	//      fp_one   or  fp_spc?
	//******************************************
	get_field(fp_one, aline, ib, &ie, pval, &nvl, &ierr);
	//read(pval(1:nvl),*, err=50) h
	sscanf(pval, "%lf", &h);
	ib = ie;
	get_field(fp_one, aline, ib, &ie, pval, &nvl, &ierr);
	double p;
	sscanf(pval, "%lf", &p);
	ib = ie;
	get_field(fp_one, aline, ib, &ie, pval, &nvl, &ierr);
	double s;
	sscanf(pval, "%lf", &s);
	ib = ie;
	get_field(fp_one, aline, ib, &ie, pval, &nvl, &ierr);

	if (nl >= MAX1D) {
		printf(" Too many layers; maximum now is %d", MAX1D);
		assert(!(nl >= MAX1D));
	}
	vp[nl][0] = p;
	if (vs1d == 1) {
		vp[nl][1] = s;
	} else {
		vp[nl][1] = p / s;
	}
//****ONE-TIME CLUDGE TO FORCE A VP/VS
//   vp(nl,2) = p/1.78
//***********************************
	z[nl] = h;
	terp[nl] = pval[0];

//   write(*,500) nl, vp(nl,1), vp(nl,2), z(nl), terp(nl)
//500    format(i4,3f10.2,2x,a1)
//   read(*,*) pause
	nl = nl + 1;
	goto a21;

//----generate the model
	a2: fp_cor = fopen(oldvfil, "wb");
	for (int n = 0; n <= 1; n++) {
		int noff = nxyzc * n;
		fprintf(fp_log, " \n");
		printf(" \n");
		fprintf(fp_log,
				"  Lay   Dep      D1      D2      V1      V2      V      ZFL     VFL\n");
		printf(
				"  Lay   Dep      D1      D2      V1      V2      V      ZFL     VFL\n");
		for (k = 0; k < nzc; k++) {
			int koff = nxyc * k + noff;
			double zg = gz[k];
			int ik;
			for (ik = 1; ik < nl; ik++) {
				if (z[ik] > zg)
					break;
			}
			ik--;
			double v;
			if (terp[ik] == 'I') {
				int zk = z[ik];
				int hz = z[ik + 1] - zk;
				double fz = (zg - zk) / hz;
				v = (1.0f - fz) * vp[ik][n] + fz * vp[ik + 1][n];
			} else {
				v = vp[ik][n];
			}
//----flatten this wavespeed if necessary
			double zfl, vfl;
			if (iflat == 1) {
				zfl = flatz(zg);
				vfl = flatvel(v, zfl);
			} else {
				zfl = zg;
				vfl = v;
			}
			fprintf(fp_log, "%4d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n",
					(k + 1), zg, z[ik], z[ik + 1], vp[ik][n], vp[ik + 1][n], v,
					zfl, vfl);
			printf("%4d%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n", (k + 1), zg,
					z[ik], z[ik + 1], vp[ik][n], vp[ik + 1][n], v, zfl, vfl);
			int j;
			for (j = 0; j < nyc; j++) {
				int joff = koff + nxc * j;
				for (i = 0; i < nxc; i++) {
					vsave[joff + i] = vfl;
				}
			}
		}
	}

	hdr_appender(hdr, nhbyte, head, type, syst, quant, flatten, hcomm);
	fwrite(hdr, sizeof(hdr[0]), nhbyte, fp_cor);
	fwrite(igridx, sizeof(igridx[0]), nxc - 1, fp_cor);
	fwrite(igridy, sizeof(igridy[0]), nyc - 1, fp_cor);
	fwrite(igridz, sizeof(igridz[0]), nzc - 1, fp_cor);
	fwrite(vsave, sizeof(vsave[0]), nxyzc2, fp_cor);

	goto a65;
	a50: printf("Error trying to read variable %s\n", mvals[i]);
	goto a65;
	a51: printf("Error trying to read filename %s\n", files[i]);

	a65: fclose(fp_spc);
	fclose(fp_log);
	fclose(fp_cor);
	fclose(fp_one);
	return 0;
}

double flatvel(double v, double z) {
// --- does earth-flattening correction for velocity, at depth z.
//     assumes z is the flat-earth depth (already corrected). 12/87 gaa.
//     a different transform is done for blocks, which have integrated
//     average velocities
	double r = 6371.00;
	return exp(z / r) * v;
}

double flatz(double z) {
// --- does earth-flattening correction from depth in spherical earth to
//     depth in flat-earth.  this preserves travel-times if the velocity
//     correction is also used.  see chapman (1973).  gaa 12/87.
	double r = 6371.00;
	return r * log(r / (r - z)) / log(exp(1));
}

double uflatz(double z) {
// --- undoes earth-flattening correction from depth in spherical earth to
//     depth in flat-earth.

	double r = 6371.00;
	return r * (1. - exp(-z / r));
}
