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
#include "common/shared_variables.h"
#include "FDtomo/make1d.h"
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

float vp[MAX1D][2], z[MAX1D];
//float vsave[nxyzcm2];

char VERSION[10] = "2004.0909\0";
double rearth = 6371.0f, degrad = 0.017453292f, hpi = 1.570796f;

char terp[MAX1D + 1];

//----header stuff
char head[5], type[5], syst[5];
char quant[5];
char flatten[5];
char hcomm[101];

//----fxs, fys, and fzs are not used in wavespeed models, so just set to zero
double fxs = 0.0, fys = 0.0, fzs = 0.0;
double axo, ayo, azo;
int nxh, nyh, nzh;
//--------------------------------------
char hdr[nhbyte + 1];

int lenhead = nhbyte;

FILE *fp_log;
FILE *fp_one;

float flatvel(float, float);
float uflatz(float);
float flatz(float);
char * dtoa(char *, double, int);

MAKE1D_DATA *make1d(SPEC spec) {
	MAKE1D_DATA *make1d_data = (MAKE1D_DATA *)malloc(sizeof(MAKE1D_DATA));
	memset(make1d_data, 0, sizeof(MAKE1D_DATA));

	//initialize variable
	int nxc = spec.grid.nxc, nyc = spec.grid.nyc, nzc = spec.grid.nzc, nx = spec.grid.nx,
	    ny = spec.grid.ny, nz = spec.grid.nz;
	double h = spec.grid.h, x0 = spec.grid.x0, *y = spec.grid.y, 
	z0 = spec.grid.z0, dq = spec.grid.dq, df = spec.grid.df, x00 = spec.grid.x00, y00 = spec.grid.y00;
	int *igridx = spec.grid.igridx, *igridy = spec.grid.igridy, *igridz = spec.grid.igridz;
	double dx = spec.grid.dx, dy = spec.grid.dy, dz = spec.grid.dz;


	float *gx = spec.grid.gx, *gy = spec.grid.gy, *gz = spec.grid.gz;

	double clat = spec.clat, clon = spec.clon, cz = spec.cz;
	float az = spec.az, azmod = spec.azmod;
	int iflat = spec.iflat, isph = spec.isph, vs1d = spec.vs1d;
	char spec_file[MAXSTRLEN];
	sscanf(spec.spec_file, "%s", spec_file);
	
	char aline[MAXSTRLEN + 1];
    int len, ierr;
	char pval[MAXSTRLEN + 1];

	int i,k;
	int ib = 0, ie = 0, lenv = 0, nvl = 0;
	
	int nxyc = nxc * nyc;
	int nxyzc = nxyc * nzc;
	int nxyzc2 = nxyzc * 2;

	
	int nxy = nx * ny;
	int nxyz = nxy * nz;

	int lengrd = 4 * (nxc + nyc + nzc - 3);
	int lenrec = lenhead + lengrd + 4 * nxyzc2;

	fp_log = fopen("make1d.log", "w");
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
		fprintf(fp_log, "  Fine Longitude spacing (df):   %.17E\n", dx);
		fprintf(fp_log, "  Fine Latidtude spacing (dq):   %.17E\n", dy);
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
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of Y coarse grid nodes: %12d\n", nyc);
	fprintf(fp_log, "  Y coarse grid node spacing: \n");
	for (int iii = 0; iii < nyc - 1; iii++) {
		fprintf(fp_log, "% 4d", igridy[iii]);
		if (iii % 10 == 9) {
			fprintf(fp_log, "\n");
		}
	}
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
	fprintf(fp_log, " One-D ASCII model file: %-40s\n", spec.onedfil);
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Output file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " New Wavespeed model file: %-40s\n", spec.oldvfil);
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
	sprintf(hcomm, "Output from make1d.c using %s", spec.onedfil);
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

	fp_one = fopen(spec.onedfil, "r");
	if(!fp_one) {
        printf("Error on opening fp_one(%s)\n", spec.onedfil);
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
	float p = 0;
	sscanf(pval, "%f", &p);
	ib = ie;
	get_field(fp_one, aline, ib, &ie, pval, &nvl, &ierr);
	float s = 0;
	sscanf(pval, "%f", &s);
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
	a2: 
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
			float zg = gz[k];
			int ik;
			for (ik = 1; ik < nl; ik++) {
				if (z[ik] > zg)
					break;
			}
			ik--;
			float v = 0;
			if (terp[ik] == 'I') {
				int zk = z[ik];
				int hz = z[ik + 1] - zk;
				float fz = (zg - zk) / hz;
				v = (1.0f - fz) * vp[ik][n] + fz * vp[ik + 1][n];
			} else {
				v = vp[ik][n];
			}
//----flatten this wavespeed if necessary
			float zfl = 0, vfl = 0;
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
					make1d_data->vsave[joff + i] = vfl;
				}
			}
		}
	}

	hdr_appender(hdr, nhbyte, head, type, syst, quant, flatten, hcomm);

	memcpy(make1d_data->hdr, hdr, strlen(hdr)+1);
	memcpy(make1d_data->igridx, igridx, 4 * (nxc - 1));
	memcpy(make1d_data->igridy, igridy, 4 * (nyc - 1));
	memcpy(make1d_data->igridz, igridz, 4 * (nzc - 1));
	//memcpy(make1d_data.vsave, vsave, sizeof(vsave));
	printf("%f\n", make1d_data->vsave[19403]);

	fclose(fp_log);
	fclose(fp_one);
	return make1d_data;
}

float flatvel(float v, float z) {
// --- does earth-flattening correction for velocity, at depth z.
//     assumes z is the flat-earth depth (already corrected). 12/87 gaa.
//     a different transform is done for blocks, which have integrated
//     average velocities
	float r = 6371.00f;
	return expf(z / r) * v;
}

float flatz(float z) {
// --- does earth-flattening correction from depth in spherical earth to
//     depth in flat-earth.  this preserves travel-times if the velocity
//     correction is also used.  see chapman (1973).  gaa 12/87.
	float r = 6371.00f;
	return r * logf(r / (r - z)) / logf(expf(1));
}

float uflatz(float z) {
// --- undoes earth-flattening correction from depth in spherical earth to
//     depth in flat-earth.

	float r = 6371.00f;
	return r * (1. - expf(-z / r));
}

int OUTPUT_MAKE1D(MAKE1D_DATA *maked1d_data, SPEC spec){
	FILE *fp_cor;
	fp_cor = fopen(spec.oldvfil, "wb");
	if (!fp_cor) {
		printf("file create error: %s\n", spec.oldvfil);
		assert(0);
	}
	
	int nxc = spec.grid.nxc, nyc = spec.grid.nyc, nzc = spec.grid.nzc;
	int nxyc = nxc * nyc;
	int nxyzc = nxyc * nzc;
	int nxyzc2 = nxyzc * 2;
	fwrite(maked1d_data->hdr, sizeof(maked1d_data->hdr[0]), nhbyte, fp_cor);
	fwrite(maked1d_data->igridx, sizeof(maked1d_data->igridx[0]), nxc - 1, fp_cor);
	fwrite(maked1d_data->igridy, sizeof(maked1d_data->igridy[0]), nyc - 1, fp_cor);
	fwrite(maked1d_data->igridz, sizeof(maked1d_data->igridz[0]), nzc - 1, fp_cor);
	fwrite(maked1d_data->vsave, sizeof(maked1d_data->vsave[0]), nxyzc2, fp_cor);
	return 0;
}

