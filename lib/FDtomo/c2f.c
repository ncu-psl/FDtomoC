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
#include "common/gridspec.h"
#include "common/parseprogs.h"
#include "common/string_process.h"
#include "common/shared_variables.h"
#include "common/read_spec.h"


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

double gx[nxcm], gy[nxcm], gz[nxcm];
float vp[nxyzcm2];

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

char hdr[nhbyte + 1];

int len_head = nhbyte;

int nxyc, nxyzc, nxyzc2;

void find_vel(double, double, double, double *, int);

int c2f(char *file_parameter) {
	char spec_file[MAXSTRLEN + 1];
	char pval[MAXSTRLEN + 1];
	char aline[MAXSTRLEN + 1];
	char varname[MAXSTRLEN + 1], parval[MAXSTRLEN + 1];
	sscanf(file_parameter, "%s", spec_file);
	spec_file[MAXSTRLEN] = '\0';

	FILE *fp_spc = fopen(spec_file, "r");
	if(!fp_spc) {
		printf("(Error in c2f.c)read fp_spc file error.\n");
		assert(0);
	}
	int len, ierr;
	int i = 0;
//c-----Grid specs

	int ib = 0, ie = 0, lenv = 0, nvl = 0;
	rewind(fp_spc);
	a11: get_line(fp_spc, aline, &ierr);
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

//c-----end of optional parameters

	nxyc = nxc * nyc;
	nxyzc = nxyc * nzc;
	nxyzc2 = nxyzc * 2;

	nx = 1, ny = 1, nz = 1;
	gx[0] = x0;
	gy[0] = y[0];
	gz[0] = z0;


// c----this is for spherical (make an option at some point)
// c        dq = h/rearth
// c        df = dabs(h/(rearth*dsin(y0)))

	for (int i = 1; i < nxc; i++) {
		nx = nx + igridx[i - 1];
		gx[i] = gx[i - 1] + h * igridx[i - 1];
	}

	for (int i = 1; i < nyc; i++) {
		ny = ny + igridy[i - 1];
		gy[i] = gy[i - 1] + h * igridy[i - 1];
	}

	for (int i = 1; i < nzc; i++) {
		nz = nz + igridz[i - 1];
		gz[i] = gz[i - 1] + h * igridz[i - 1];
	}
//	c----dimension check
	if (nx > nxm) {
		printf(" nx(%d) is too large, maximum is: %d\n",nx , nxm);
		assert(0);
	}
	if (ny > nym) {
		printf(" ny is too large, maximum is: %d\n", nym);
		assert(0);
	}
	if (nz > nzm) {
		printf(" nz is too large, maximum is: %d\n", nzm);
		assert(0);
	}

	int nxy = nx * ny;
	int nxyz = nxy * nz;
	int nxyz2 = nxyz * 2;

	int lenlog;
	char logfile[80 + 1];
	sprintf(logfile, "c2f.log%d", ittnum);
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
			spec_file);
	fprintf(fp_log, "  \n");
	fprintf(fp_log, "  Iteration counter:%19d\n", ittnum);
	fprintf(fp_log, "  \n");
	char out_str[MAXSTRLEN];
	dtoa(out_str, x0, 18);
	fprintf(fp_log, "  Cartesian X origin (x0):   %s     \n", out_str);
	dtoa(out_str, y[0], 18);
	fprintf(fp_log, "  Cartesian Y origin (y0):   %s     \n", out_str);
	dtoa(out_str, z0, 18);
	fprintf(fp_log, "  Cartesian Z origin (z0):   %s     \n", out_str);
	dtoa(out_str, h, 18);
	fprintf(fp_log, "  Fine grid spacing:   %s     \n", out_str);
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of X coarse grid nodes: %12d\n", nxc);
	fprintf(fp_log, "  X coarse grid node spacing: \n");
	//write(lout, "(10i4)")(igridx(i), i = 1, nxc - 1)
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
//write(lout,"(10i4)") (igridy(i),i=1,nyc-1)
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
//write(lout,"(10i4)") (igridz(i),i=1,nzc-1)
	for (int iii = 0; iii < nzc - 1; iii++) {
		fprintf(fp_log, "% 4d", igridz[iii]);
		if (iii % 10 == 9) {
			fprintf(fp_log, "\n");
		}
	}
	fprintf(fp_log, "\n");
	fprintf(fp_log, " \n");
	//fprintf(fp_log, "  Number of X fine grid nodes: %12d\n", nx+1);
	//fprintf(fp_log, "  Number of Y fine grid nodes: %12d\n", ny+1);
	//fprintf(fp_log, "  Number of Z fine grid nodes: %12d\n", nz+1);
	fprintf(fp_log, "  Number of X fine grid nodes: %12d\n", nx);
	fprintf(fp_log, "  Number of Y fine grid nodes: %12d\n", ny);
	fprintf(fp_log, "  Number of Z fine grid nodes: %12d\n", nz);
	fprintf(fp_log, " \n");
	//fprintf(fp_log, "  Total Number of fine grid nodes:%15d\n", (nx+1)*(ny+1)*(nz+1));
	fprintf(fp_log, "  Total Number of fine grid nodes:%15d\n", nx * ny * nz);
	fprintf(fp_log, "  Total Number of coarse grid nodes:%13d\n", nxyzc);
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Input file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Current Wavespeed model file: %-40s\n", oldvfil);
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Output file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Grid output file:             %-40s\n", tgrdfil);
	fprintf(fp_log, " P&S fine model:               %-40s\n", finevel);
	fprintf(fp_log, " P only fine model (*.pvel):   %-40s\n", finevel);
	fprintf(fp_log, " S only fine model (*.svel):   %-40s\n", finevel);
	fprintf(fp_log, " \n");
//c---temporary output of grid for testing.  This will not be coorect
//c   for spherical coordinates.
	FILE *fp_tgd = fopen(tgrdfil, "w");
	if (!fp_tgd) {
		printf("(Error in c2f.c)create file error.\n");
		assert(0);
	}
	for (i = 0; i < nxc; i++) {
		dtoa(out_str, gx[i], 10);
		fprintf(fp_tgd, "  %s    ", out_str);
		dtoa(out_str, gy[0], 10);
		fprintf(fp_tgd, "  %s    \n", out_str);
		dtoa(out_str, gx[i], 10);
		fprintf(fp_tgd, "  %s    ", out_str);
		dtoa(out_str, gy[nyc - 1], 10);
		fprintf(fp_tgd, "  %s    \n", out_str);
	}

	for (i = 0; i < nyc; i++) {
		dtoa(out_str, gx[0], 10);
		fprintf(fp_tgd, "  %s    ", out_str);
		dtoa(out_str, gy[i], 10);
		fprintf(fp_tgd, "  %s    \n", out_str);
		dtoa(out_str, gx[nxc - 1], 10);
		fprintf(fp_tgd, "  %s    ", out_str);
		dtoa(out_str, gy[i], 10);
		fprintf(fp_tgd, "  %s    \n", out_str);
	}

// c---tempoary output for plotting
// c	do j = 1, nyc
// c	  do i = 1, nxc
// c	    write(23,1234) i, j,  gx(i)*1000.0, gy(j)*1000.0
// c1234	    format(2i4, ' Q ', 2f12.2)
// c	  enddo
// c	enddo

	char *cfile = oldvfil;
	FILE *fp_cor = fopen(cfile, "rb");
	if (!fp_cor) {
		printf("(Error in c2f.c)read file error.\n");
		assert(0);
	}
	fread(hdr, sizeof(hdr[0]), nhbyte, fp_cor);
	
	char *offset = hdr;
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
	sscanf(offset, "%100s", hcomm);
	
	for (int i = 0; i < nxyzm2; i++)
		vsave[i] = 0.;

//c---verify that this is a valid header
	int len_rec;
	if (strcmp(head, "HEAD") != 0) {
		printf(
				"File does not contain valid header...attempting headerless read \n");
		fprintf(fp_log,
				"File does not contain valid header...attempting headerless read \n");
		len_rec = 4 * nxyzc2;
		rewind(fp_cor);
		if (!fp_cor) {
			printf("(Error in c2f.c)read file error.\n");
			assert(0);
		}
		fread(vp, sizeof(vp[0]), nxyzc2, fp_cor);
//c---set trial header values
		strcpy(head, "HEAD");
		strcpy(syst, "CART");
		strcpy(quant, "BMOD");
		strcpy(flatten, "NOFL");
		clath = clat;
		clonh = clon;
		czh = cz;
		azh = azmod;
		axo = x0;
		ayo = y[0];
		azo = z0;
		dxh = h;
		dyh = h;
		dzh = h;
		nxh = nxc;
		nyh = nyc;
		nzh = nzc;
	} else {
		if (strcmp(type, "CORS") != 0) {
			printf(" WARNING: input mesh does not appear to be COARSE: %s\n",
					type);
			fprintf(fp_log,
					" WARNING: input mesh does not appear to be COARSE: %s\n",
					type);
		}
		if (strcmp(quant, "BMOD") != 0) {
			printf(" WARNING: file does not appear to be a valid type: %s\n",
					quant);
			fprintf(fp_log,
					" WARNING: file does not appear to be a valid type: %s\n",
					quant);
		}

		printf(" Reading in Coarse Mesh ...");
		int len_grd = 4 * (nxc + nyc + nzc - 3);
		len_rec = len_head + len_grd + 4 * nxyzc2;

		fread(igridx, sizeof(igridx[0]), nxc - 1, fp_cor);
		fread(igridy, sizeof(igridy[0]), nyc - 1, fp_cor);
		fread(igridz, sizeof(igridz[0]), nzc - 1, fp_cor);
		fread(vp, sizeof(vp[0]), nxyzc2, fp_cor);
	}
	printf(".. Done \n");

// c---tempoary output for plotting
// c   i1 = 0
// c   do k = 1, nzc
// c     do j = 1, nyc
// c       do i = 1, nxc
// c       i1 = i1 + 1
// c       write(24,1244) i, j, gz(k)*1000.0,  vp(i1), vp(i1+nxyzc),
// c     +       vp(i1)/vp(i1+nxyzc), gx(i)*1000.0, gy(j)*1000.0
// c1244       format(2i4, f12.2, 3f6.2, 8f12.2)
// c       enddo
// c     enddo
// c   enddo

// c****ONE TIME CLUDGE TO FORCE CONSTANT VPVS
// c   do i = 1, nxyzc
// c     vp(i+nxyzc) = vp(i)/1.78
// c   enddo
// c*****END****

// c----convert to slowness
	for (i = 0; i < nxyzc2; i++) {
		vp[i] = 1. / vp[i];
	}

	printf("Interpolating ...");
// c----interpolation
	int iph, joff, nk2, ij = 0;
	double z, yy, x, v;
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
					find_vel(x, yy, z, &v, iph);
//c-----convert back to velocity
					vsave[ij] = 1. / v;
					ij++;
				}
			}
		}
	}
	printf(".. Done \n");

//c---redefine mesh type to be "FINE"
	strcpy(type, "FINE");
	for (i = 0; i < 100; i++)
		hcomm[i] = ' ';

	sprintf(hcomm, "Output from c2f.c Version %s using %s", VERSION, oldvfil);
	nxh = nx;
	nyh = ny;
	nzh = nz;
//c---write out

	char ffile[MAXSTRLEN + 1];
	strcpy(ffile, finevel);
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
	hdr_appender(hdr, nhbyte, head, type, syst, "BMOD", flatten, hcomm);

	fwrite(hdr, sizeof(hdr[0]), nhbyte, fp_fnw);
	fwrite(vsave, sizeof(vsave[0]), nxyz2, fp_fnw);

	FILE *fp_fnp = fopen(vpfile, "wb");
	if (!fp_fnp) {
		printf("(Error in c2f.c)write file error.\n");
		assert(0);
	}
	printf("Writing out file 2... \n");
	hdr_appender(hdr, nhbyte, head, type, syst, "VPMD", flatten, hcomm);

	fwrite(hdr, sizeof(hdr[0]), nhbyte, fp_fnp);
	fwrite(vsave, sizeof(vsave[0]), nxyz, fp_fnp);

	FILE *fp_fns = fopen(vsfile, "wb");
	if (!fp_fns) {
		printf("(Error in c2f.c)write file error.\n");
		assert(0);
	}
	printf("Writing out file 3... \n");
	hdr_appender(hdr, nhbyte, head, type, syst, "VSMD", flatten, hcomm);

	fwrite(hdr, sizeof(hdr[0]), nhbyte, fp_fns);
	fwrite(vsave + nxyz, sizeof(vsave[0]), nxyz2 - nxyz, fp_fns);

	a65: 
	fclose(fp_spc);
	fclose(fp_log);
	fclose(fp_cor);
	fclose(fp_tgd);
	fclose(fp_fnw);
	fclose(fp_fnp);
	fclose(fp_fns);
	return 0;
}

void find_vel(double x, double y, double z, double *v, int iph) {
// c-----find a velocity at x, y, z
	int i;
	for (i = 1; i < nxc; i++) {
		if (gx[i] > x)
			break;
	}
	if (i == nxc)
		i--;
	i--;
	double xi = gx[i];
	double hx = gx[i + 1] - xi;

	int j;
	for (j = 1; j < nyc; j++) {
		if (gy[j] > y)
			break;
	}
	if (j == nyc)
		j--;
	j--;
	double yj = gy[j];
	double hy = gy[j + 1] - yj;

	int k;
	for (k = 1; k < nzc; k++) {
		if (gz[k] > z)
			break;
	}
	if (k == nzc)
		k--;
	k--;
	double zk = gz[k];
	double hz = gz[k + 1] - zk;

	int nk = nxc * nyc * k;
	int nj = nxc * j;
	int nk2 = nxc * nyc * (k + 1);
	int nj2 = nxc * (j + 1);
	double fx = (x - xi) / hx;
	double fy = (y - yj) / hy;
	double fz = (z - zk) / hz;
	i = i + iph * nxyzc;
	double tv = 0.;
	tv = (1. - fx) * (1. - fy) * (1. - fz) * vp[nk + nj + i];
	tv = tv + fx * (1. - fy) * (1. - fz) * vp[nk + nj + i + 1];
	tv = tv + (1. - fx) * fy * (1. - fz) * vp[nk + nj2 + i];
	tv = tv + (1. - fx) * (1. - fy) * fz * vp[nk2 + nj + i];
	tv = tv + fx * fy * (1. - fz) * vp[nk + nj2 + i + 1];
	tv = tv + fx * (1. - fy) * fz * vp[nk2 + nj + i + 1];
	tv = tv + (1. - fx) * fy * fz * vp[nk2 + nj2 + i];
	tv = tv + fx * fy * fz * vp[nk2 + nj2 + i + 1];
	*v = tv;
	return;
}
