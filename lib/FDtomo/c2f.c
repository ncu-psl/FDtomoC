// c
// c	Version 2005.0809
// c
// c	Program to convert from coarse grid velocity files to fine grid
// c	velocity files.  Useful for input to punch and other programs requring
// c	a fine grid.
// c
// c	Note that in this version we interpolate slownesses rather than velocities
// c	since the eikonal sofware (punch.c) and related sofware works in slowness.
// c
// c	This version uses headers.
// c	
// c	Input files
// c
// c       fn      Name   Uname   U#  Comments
// c        5     oldvfil luncor  12  Coarse grid wavespeed model file
// c
// c	Output files
// c
// c       fn      Name   		Uname   U#  Comments
// c	26	tgrdfil		luntgd	34  Grid output for plotting and test purposes
// c	 6	finevel		lunfnw	13  ffile   Fine grid (P&S)
// c	 6	finevel.pvel	lunfnp	14  vpfile  Fine grid (P only)
// c	 6	finevel.svel	lunfns	15  vsfile  Fine grid (S only)
// c
// c	Note on adaptation to spherical system:
// c	Because we are interpolating over independant axes (x, y, z or f, q, r) the
// c	dimensions we choose are not important.  Hence we can use h for all directions
// c	rather than worrying about dq and df.   The only difference comes in computing
// c	the actual grid locations for purposes of mapping.  The grid output for spherical
// c	systems is not correct in this version of c2f.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <limits.h>

#include "common/environment_setting.h"
#include "common/parameter.h"
#include "common/parseprogs.h"
#include "common/string_process.h"
#include "common/gridspec.h"
#include "FDtomo/c2f.h"

#define MAX1D 1000
#define MAXSTRLEN 132
// c--mustv is the number of variables that must be assigned in the parameter
// c       file in order for this program to run correctly.  The names of the
// c       variables are specified in the mvals string array below.
// c--mustf is the number of files that must be attached; see the files array below.
#define MUSTV  4
#define MUSTF  3

// c---number of 4 byte words in the header
#define nhbyte 58 * 4
float *vp;

float vsave[nxyzm2];

char cfile[MAXSTRLEN + 1], ffile[MAXSTRLEN + 1], vpfile[MAXSTRLEN + 1], vsfile[MAXSTRLEN + 1];

//c----header stuff
char head[5], type[5], syst[5];
char quant[5];
char flatten[5];
char hcomm[101];

double fxs, fys, fzs;
double clath, clonh, czh;
double axo, ayo, azo, dxh, dyh, dzh;
float azh;
int nxh, nyh, nzh;

char hdr[120];

int len_head = nhbyte;

int nxyc, nxyzc, nxyzc2;

void find_vel(float, float, float, float *, int, GRID grid);

C2F_DATA *c2f(SPEC spec, MAKE1D_DATA *MAKE1D) {
	C2F_DATA *c2f_data = (C2F_DATA *)malloc(sizeof(C2F_DATA));
	c2f_data->vpfile = (VELFILE *)malloc(sizeof(VELFILE));
	c2f_data->vsfile = (VELFILE *)malloc(sizeof(VELFILE));



	//initialize variable
	int nxc = spec.grid.nxc, nyc = spec.grid.nyc, nzc = spec.grid.nzc, nx = spec.grid.nx,
	    ny = spec.grid.ny, nz = spec.grid.nz;
	double h = spec.grid.h, x0 = spec.grid.x0, *y = spec.grid.y, 
	z0 = spec.grid.z0, dq = spec.grid.dq, df = spec.grid.df, x00 = spec.grid.x00, y00 = spec.grid.y00;
	int *igridx = spec.grid.igridx, *igridy = spec.grid.igridy, *igridz = spec.grid.igridz;

	float *gx = spec.grid.gx, *gy = spec.grid.gy, *gz = spec.grid.gz;

	double clat = spec.clat, clon = spec.clon, cz = spec.cz;
	float az = spec.az, azmod = spec.azmod;
	int iflat = spec.iflat, isph = spec.isph, vs1d = spec.vs1d;
	int ittnum = spec.ittnum;
	double rearth = 6371.0, degrad = 0.017453292, hpi = 1.570796;

	char spec_file[MAXSTRLEN];
	sscanf(spec.spec_file, "%s", spec_file);
	char VERSION[10] = "2005.0809\0";

	int len, ierr;
	int i = 0, k;

	nxyc = nxc * nyc;
	nxyzc = nxyc * nzc;
	nxyzc2 = nxyzc * 2;
	
	int nxy = nx * ny;
	int nxyz = nxy * nz;
	int nxyz2 = nxyz * 2;

	int lenlog;
	
	for (int i = 0; i < nxyzm2; i++)
		vsave[i] = 0.;

//c---verify that this is a valid header
	struct vhead *headin = &MAKE1D->head;
	int len_rec;

	int len_grd = 4 * (nxc + nyc + nzc - 3);
	len_rec = len_head + len_grd + 4 * nxyzc2;
	igridx = MAKE1D->igridx;
	igridy = MAKE1D->igridy;
	igridz = MAKE1D->igridz;
	vp = MAKE1D->vsave;
	printf(".. Done \n");

// c----convert to slowness
	for (i = 0; i < nxyzc2; i++) {
		vp[i] = 1. / vp[i];
	}

	printf("Interpolating ...");
// c----interpolation
	int iph, joff, nk2, ij = -1;
	float z, yy, x, v;
	for (int n = 0; n <= 1; n++) {
		iph = n;
		joff = ny * nz * n;
		for (k = 0; k < nz; k++) {
			z = z0 + k * h;
			nk2 = k * ny + joff;
			for (int j = 0; j < ny; j++) {
				yy = y[0] + j * h;
				for (i = 0; i < nx; i++) {
					x = x0 + i * h;
					find_vel(x, yy, z, &v, iph, spec.grid);
//c-----convert back to velocity
					ij++;
					vsave[ij] = 1.f / v;
				}
			}
		}
	}
	printf(".. Done \n");

//c---write out

	char ffile[MAXSTRLEN + 1];
	strcpy(ffile, spec.finevel);
	int leng_thp = strlen(ffile);
	for (int i = 0; i < strlen(ffile); i++) {
		if (ffile[i] == ' ') {
			leng_thp = i;
			break;
		}
	}
	ffile[leng_thp] = '\0';

	sprintf(vpfile, "%s.pvel", ffile);
	sprintf(vsfile, "%s.svel", ffile);


	FILE *fp_fnw = fopen(ffile, "wb");
	if (!fp_fnw) {
		printf("(Error in c2f.c)write file error.\n");
		assert(0);
	}
	printf("Writing out file 1... \n");

	fwrite(vsave, sizeof(vsave[0]), nxyz2, fp_fnw);

	memcpy(c2f_data->vpfile->filename, vpfile, strlen(vpfile)+1);
	memcpy(c2f_data->vpfile->vsave, vsave, sizeof(vsave[0])*nxyz);

	memcpy(c2f_data->vsfile->filename, vsfile, strlen(vsfile)+1);
	memcpy(c2f_data->vsfile->vsave, vsave + nxyz, sizeof(vsave[0]) * (nxyz2 - nxyz));

	a65: 
	fclose(fp_fnw);
	return c2f_data;
}

void find_vel(float x, float y, float z, float *v, int iph, GRID grid) {
// c-----find a velocity at x, y, z
	int i;
	int nxc = grid.nxc, nyc = grid.nyc, nzc = grid.nzc;
	float *gx = grid.gx, *gy = grid.gy, *gz = grid.gz;

	for (i = 1; i < nxc; i++) {
		if (gx[i] > x)
			break;
	}
	if (i == nxc)
		i--;
	i--;

	float xi = gx[i];
	float hx = gx[i + 1] - xi;
	int j;
	for (j = 1; j < nyc; j++) {
		if (gy[j] > y)
			break;
	}
	if (j == nyc)
		j--;
	j--;

	float yj = gy[j];
	float hy = gy[j + 1] - yj;
	int k;
	for (k = 1; k < nzc; k++) {
		if (gz[k] > z)
			break;
	}
	if (k == nzc)
		k--;
	k--;
	float zk = gz[k];
	float hz = gz[k + 1] - zk;

	int nk = nxc * nyc * k;
	int nj = nxc * j;
	int nk2 = nxc * nyc * (k + 1);
	int nj2 = nxc * (j + 1);
	float fx = (x - xi) / hx;
	float fy = (y - yj) / hy;
	float fz = (z - zk) / hz;
	i = i + iph * nxyzc;
	float tv = 0.;
	tv = (1.f - fx) * (1.f - fy) * (1.f - fz) * vp[nk + nj + i];
	tv = tv + fx * (1.f - fy) * (1.f - fz) * vp[nk + nj + i + 1];
	tv = tv + (1.f - fx) * fy * (1.f - fz) * vp[nk + nj2 + i];
	tv = tv + (1.f - fx) * (1.f - fy) * fz * vp[nk2 + nj + i];
	tv = tv + fx * fy * (1.f - fz) * vp[nk + nj2 + i + 1];
	tv = tv + fx * (1.f - fy) * fz * vp[nk2 + nj + i + 1];
	tv = tv + (1.f - fx) * fy * fz * vp[nk2 + nj2 + i];
	tv = tv + fx * fy * fz * vp[nk2 + nj2 + i + 1];
	*v = tv;
	return;
}

OUTPUT_C2F(C2F_DATA *c2f_data, SPEC spec){
	int nxc = spec.grid.nxc, nyc = spec.grid.nyc, nzc = spec.grid.nzc, nx = spec.grid.nx,
	    ny = spec.grid.ny, nz = spec.grid.nz;

	int nxy = nx * ny;
	int nxyz = nxy * nz;
	int nxyz2 = nxyz * 2;

	FILE *fp_fnp = fopen(c2f_data->vpfile->filename, "wb");
	if (!fp_fnp) {
		printf("(Error in c2f.c)write file error.\n");
		assert(0);
	}
	printf("Writing out file 2... \n");
	fwrite(&c2f_data->vpfile->head, sizeof(char), nhbyte, fp_fnp);
	fwrite(c2f_data->vpfile->vsave, sizeof(c2f_data->vpfile->vsave[0]), nxyz, fp_fnp);
	
	FILE *fp_fns = fopen(c2f_data->vsfile->filename, "wb");
	if (!fp_fns) {
		printf("(Error in c2f.c)write file error.\n");
		assert(0);
	}
	printf("Writing out file 3... \n");
	fwrite(&c2f_data->vsfile->head, sizeof(char), nhbyte, fp_fns);
	fwrite(c2f_data->vsfile->vsave, sizeof(c2f_data->vsfile->vsave[0]), nxyz2 - nxyz, fp_fns);


}

int LOG_C2F(SPEC spec){
	char logfile[80 + 1];
	sprintf(logfile, "c2f.log%d", spec.ittnum);
	FILE *fp_log = fopen(logfile, "w");
	if (!fp_log) {
		printf("create %s file error.\n", logfile);
		assert(0);
	}

	fprintf(fp_log, "  \n");
	fprintf(fp_log,
			" *************************************************************** \n");

	fprintf(fp_log, "          Parameters Set For This Run of c2f.c\n");

	fprintf(fp_log, "  \n");
	fprintf(fp_log, " c2f VERSION: %s\n", VERSION);
	fprintf(fp_log, "  \n");
	fprintf(fp_log, " Current parameter specification file: %-40s\n",
			spec.spec_file);
	fprintf(fp_log, "  \n");
	fprintf(fp_log, "  Iteration counter:%19d\n", spec.ittnum);
	fprintf(fp_log, "  \n");
	char out_str[MAXSTRLEN];
	dtoa(out_str, spec.grid.x0, 18);
	fprintf(fp_log, "  Cartesian X origin (x0):   %s     \n", out_str);
	dtoa(out_str, spec.grid.y[0], 18);
	fprintf(fp_log, "  Cartesian Y origin (y0):   %s     \n", out_str);
	dtoa(out_str, spec.grid.z0, 18);
	fprintf(fp_log, "  Cartesian Z origin (z0):   %s     \n", out_str);
	dtoa(out_str, spec.grid.h, 18);
	fprintf(fp_log, "  Fine grid spacing:   %s     \n", out_str);
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of X coarse grid nodes: %12d\n", spec.grid.nxc);
	fprintf(fp_log, "  X coarse grid node spacing: \n");
	//write(lout, "(10i4)")(igridx(i), i = 1, nxc - 1)
	for (int iii = 0; iii < spec.grid.nxc - 1; iii++) {
		fprintf(fp_log, "% 4d", spec.grid.igridx[iii]);
		if (iii % 10 == 9) {
			fprintf(fp_log, "\n");
		}
	}
	fprintf(fp_log, "\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of Y coarse grid nodes: %12d\n", spec.grid.nyc);
	fprintf(fp_log, "  Y coarse grid node spacing: \n");
//write(lout,"(10i4)") (igridy(i),i=1,nyc-1)
	for (int iii = 0; iii < spec.grid.nyc - 1; iii++) {
		fprintf(fp_log, "% 4d", spec.grid.igridy[iii]);
		if (iii % 10 == 9) {
			fprintf(fp_log, "\n");
		}
	}
	fprintf(fp_log, "\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of Z coarse grid nodes: %12d\n", spec.grid.nzc);
	fprintf(fp_log, "  Z coarse grid node spacing: \n");
//write(lout,"(10i4)") (igridz(i),i=1,nzc-1)
	for (int iii = 0; iii < spec.grid.nzc - 1; iii++) {
		fprintf(fp_log, "% 4d", spec.grid.igridz[iii]);
		if (iii % 10 == 9) {
			fprintf(fp_log, "\n");
		}
	}
	fprintf(fp_log, "\n");
	fprintf(fp_log, " \n");
	//fprintf(fp_log, "  Number of X fine grid nodes: %12d\n", nx+1);
	//fprintf(fp_log, "  Number of Y fine grid nodes: %12d\n", ny+1);
	//fprintf(fp_log, "  Number of Z fine grid nodes: %12d\n", nz+1);
	fprintf(fp_log, "  Number of X fine grid nodes: %12d\n", spec.grid.nx);
	fprintf(fp_log, "  Number of Y fine grid nodes: %12d\n", spec.grid.ny);
	fprintf(fp_log, "  Number of Z fine grid nodes: %12d\n", spec.grid.nz);
	fprintf(fp_log, " \n");
	//fprintf(fp_log, "  Total Number of fine grid nodes:%15d\n", (nx+1)*(ny+1)*(nz+1));
	fprintf(fp_log, "  Total Number of fine grid nodes:%15d\n", spec.grid.nx * spec.grid.ny * spec.grid.nz);
	fprintf(fp_log, "  Total Number of coarse grid nodes:%13d\n", nxyzc);
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Input file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Current Wavespeed model file: %-40s\n", spec.oldvfil);
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Output file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Grid output file:             %-40s\n", spec.tgrdfil);
	fprintf(fp_log, " P&S fine model:               %-40s\n", spec.finevel);
	fprintf(fp_log, " P only fine model (*.pvel):   %-40s\n", spec.finevel);
	fprintf(fp_log, " S only fine model (*.svel):   %-40s\n", spec.finevel);
	fprintf(fp_log, " \n");
//c---temporary output of grid for testing.  This will not be coorect
//c   for spherical coordinates.
	FILE *fp_tgd = fopen(spec.tgrdfil, "w");
	if (!fp_tgd) {
		printf("(Error in c2f.c)create file error.\n");
		assert(0);
	}
	for (int i = 0; i < spec.grid.nxc; i++) {
		dtoa(out_str, spec.grid.gx[i], 10);
		fprintf(fp_tgd, "  %s    ", out_str);
		dtoa(out_str, spec.grid.gy[0], 10);
		fprintf(fp_tgd, "  %s    \n", out_str);
		dtoa(out_str, spec.grid.gx[i], 10);
		fprintf(fp_tgd, "  %s    ", out_str);
		dtoa(out_str, spec.grid.gy[spec.grid.nyc - 1], 10);
		fprintf(fp_tgd, "  %s    \n", out_str);
	}

	for (int i = 0; i < spec.grid.nyc; i++) {
		dtoa(out_str, spec.grid.gx[0], 10);
		fprintf(fp_tgd, "  %s    ", out_str);
		dtoa(out_str, spec.grid.gy[i], 10);
		fprintf(fp_tgd, "  %s    \n", out_str);
		dtoa(out_str, spec.grid.gx[spec.grid.nxc - 1], 10);
		fprintf(fp_tgd, "  %s    ", out_str);
		dtoa(out_str, spec.grid.gy[i], 10);
		fprintf(fp_tgd, "  %s    \n", out_str);
	}

	fclose(fp_log);
	fclose(fp_tgd);
}