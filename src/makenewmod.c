#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <limits.h>
#include <memory.h>

#include "common/parameter.h"
#include "common/gridspec.h"
#include "common/parseprogs.h"
#include "common/string_process.h"
#include "vhead.h"

#define MAX1D 1000
#define MAXSTRLEN 132

// c--mustv is the number of variables that must be assigned in the parameter
// c       file in order for this program to run correctly.  The names of the
// c       variables are specified in the mvals string array below.
// c--mustf is the number of files that must be attached; see the files array below.
#define MUSTV  3
#define MUSTF  3

//---number of 4 byte words in the header
#define nhbyte 58 * 4
float rearth = 6371.0, degrad = 0.017453292, hpi = 1.570796;
#define VERSION "2018.0625"

// c---du holds all the perturbations.  
// c   ds holds the smoothed perturbations
// c   vo holds the original model
// c   vn holds the new model
float du[nxyzcm2], ds[nxyzcm2], vo[nxyzcm2], vn[nxyzcm2];

// c---vd are the pertubations from runlsqr
float vd[maxmbl];
int jndx[maxmbl];

int main(void) {
	char spec_file[MAXSTRLEN + 1];
	char pval[MAXSTRLEN + 1];
	char aline[MAXSTRLEN + 1], bline[maxlst][MAXSTRLEN + 1];
	char varname[MAXSTRLEN + 1], parval[MAXSTRLEN + 1];
	char *mvals[MUSTV] = { "nxc\0", "nyc\0", "nzc\0" };
	char *files[MUSTF] = { "nmodfil\0", "oldvfil\0", "fmodfil\0" };
	char nmodfil[MAXSTRLEN + 1], oldvfil[MAXSTRLEN + 1], fmodfil[MAXSTRLEN + 1];
	double fxs = 6.0134700169990685e-154, fys = 6.0134700169990685e-154, fzs = 6.0134700169990685e-154;
	float azh;
	double clath, clonh, czh;
	double axo, ayo, azo, dx, dy, dz;
	int nxh, nyh, nzh;
	printf("Enter parameter specification file: ");
	scanf("%s", spec_file);
	
	spec_file[MAXSTRLEN] = '\0';
	FILE *fp_spc = fopen(spec_file, "r");
	if (!fp_spc) {
		printf("error on opening spec-file (%s)\n", spec_file);
		assert(0);
	}

// c---recover the variables needed to run this program
// c
// c       nxc, nyc, nzc      coarse dimensions of the fine mesh used in the trt tables
// c
	int len, ierr;
	for (int i = 0; i < MUSTV; i++) {
		get_vars(fp_spc, mvals[i], pval, &len, &ierr);
		if (ierr == 1) {
			printf("Error trying to read variable %s\n", mvals[i]);
			assert(0);
		}
		if (i == 0) {
			sscanf(pval, "%d", &nxc);
		} else if (i == 1) {
			sscanf(pval, "%d", &nyc);
		} else if (i == 2) {
			sscanf(pval, "%d", &nzc);
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

	for (int i = 0; i < MUSTF; i++) {
		get_vars(fp_spc, files[i], pval, &len, &ierr);
		if (ierr == 1) {
			printf("Error trying to read filename %s\n", files[i]);
			assert(0);
		}
		if (i == 0) {
			sscanf(pval, "%s", nmodfil);
			//lenf1 = len;
		} else if (i == 1) {
			sscanf(pval, "%s", oldvfil);
			//lenf2 = len;
		} else if (i == 2) {
			sscanf(pval, "%s", fmodfil);
			//lenf3 = len;
		}
	}
	// c-----------------------------------------------------------------------------
// c
// c    Default setting for some variables
// c
// c---limitu = 0 means don't apply maximum limits on perturbation. = 1 means to apply

	int istacor = 0, limitu = 0, ivpvs = 1, mavx = 3, mavy = 3, mavz = 3;
	int nsmooth = 2, ipscflg = 0, ido1d = 0, ittnum = 1;

	float dvperc = 0.05, pertscl = 1.0;
	// c---End of default values

	char stafile[MAXSTRLEN + 1], nstafil[MAXSTRLEN + 1];

	// c--Optionally read in some variables
	get_vars(fp_spc, "istacor ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &istacor);
	if (istacor != 0 && istacor != 1) {
		istacor = 0;
	}
	get_vars(fp_spc, "limitu ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &limitu);
	if (limitu != 0 && limitu != 1) {
		limitu = 0;
	}
	get_vars(fp_spc, "ivpvs ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &ivpvs);
	if (ivpvs != 0 && ivpvs != 1) {
		ivpvs = 0;
	}
	get_vars(fp_spc, "mavx ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &mavx);
	get_vars(fp_spc, "mavy ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &mavy);
	get_vars(fp_spc, "mavz ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &mavz);
	get_vars(fp_spc, "nsmooth ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &nsmooth);
	get_vars(fp_spc, "dvperc ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &dvperc);
	get_vars(fp_spc, "ipscflg ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &ipscflg);
	get_vars(fp_spc, "pertscl ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &pertscl);
	get_vars(fp_spc, "do1d ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &ido1d);
	get_vars(fp_spc, "ittnum ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &ittnum);

	// Optional files
	get_vars(fp_spc, "stafile ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", stafile);
	get_vars(fp_spc, "nstafil ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", nstafil);
	// end of optional parameters
	printf("  nxc, nyc, nzc: %12d%12d%12d\n", nxc, nyc, nzc);

	int nxyc = nxc * nyc;
	int nxyzc = nxyc * nzc;
	int nxyzc2 = nxyzc * 2;

	char logfile[80 + 1];
	sprintf(logfile, "sphfdloc.log%d", ittnum);
	FILE *fp_log = fopen(logfile, "w");
	if (!fp_log) {
		printf("create %s file error.\n", logfile);
		assert(0);
	}

	fprintf(fp_log, "  \n");
	fprintf(fp_log,
			" *************************************************************** \n");
	fprintf(fp_log, "          Parameters Set For This Run of Makenewmod.c\n");
	fprintf(fp_log, "  \n");
	fprintf(fp_log, "Sphfdloc VERSION: %s\n", VERSION);
	fprintf(fp_log, "  \n");
	fprintf(fp_log, " Current parameter specification file: %-40s\n",
			spec_file);
	fprintf(fp_log, "  \n");
	{
		char tmp[MAXSTRLEN];
		dtoa(tmp, ittnum, 18);
		fprintf(fp_log, " Iteration counter:          %s     \n", tmp);
	}
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of X fine grid nodes: %12d\n", nxc);
	fprintf(fp_log, "  Number of Y fine grid nodes: %12d\n", nyc);
	fprintf(fp_log, "  Number of Z fine grid nodes: %12d\n", nzc);
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of X fine grid nodes: %12d\n", nxc);
	fprintf(fp_log, "  Number of Y fine grid nodes: %12d\n", nyc);
	fprintf(fp_log, "  Number of Z fine grid nodes: %12d\n", nzc);
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Moving Window Length in X : %12d\n", mavx);
	fprintf(fp_log, "  Moving Window Length in Y : %12d\n", mavy);
	fprintf(fp_log, "  Moving Window Length in Z : %12d\n", mavz);
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of smoothing passes : %12d\n", nsmooth);
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  ivpvs: %12d\n", ivpvs);
	fprintf(fp_log, " \n");

	if (ipscflg == 1) {
		fprintf(fp_log, " Perturbations scaled at the outset by : %f", pertscl);
	} else {
		fprintf(fp_log, " Perturbations not scaled at the outset (pscflg = 0)");
	}
	if (limitu == 1) {
		fprintf(fp_log, " Maximum Wavespeed Percentage change allowed : %f",
				dvperc);
	} else {
		fprintf(fp_log, " No perturbation limits applied (limitu = 0) ");
	}
	if (ido1d == 1) {
		fprintf(fp_log, " Corrections assumed for 1D model (do1d = 1) ");
	} else {
		fprintf(fp_log, " Corrections assumed for 3D model (do1d = 0) ");
	}

	fprintf(fp_log, " Input file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Arr. Time Wavespeed model: %s\n", oldvfil);
	fprintf(fp_log, " Perturbation file: %s\n", nmodfil);
	fprintf(fp_log, " Existing Station List: %s\n", stafile);
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Output file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " New Wavespeed Model: %s\n", fmodfil);
	fprintf(fp_log, " New Station List : %s\n", nstafil);
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Total Number of coarse grid nodes:%13d\n", nxyzc);
	fprintf(fp_log, " \n");
	fprintf(fp_log,
			" *************************************************************** \n");
	int ima = mavx;
	int jma = mavy;
	int kma = mavz;
	if ((ima % 2) != 1 || (jma % 2) != 1 || kma % 2 != 1) {
		printf("*** ERROR: invalid size for moving average filter\n");
		assert(0);
	}
	ima /= 2;
	jma /= 2;
	kma /= 2;

// c-----open the station file; terminate with blank line or EOF
	int nstr = 0;
	float pcor[maxlst], scor[maxlst];
	memset(pcor, 0, sizeof(pcor));
	memset(scor, 0, sizeof(scor));
	if (istacor == 1) {
		FILE * fp_sta = fopen(stafile, "r");
		if (!fp_sta) {
			printf("file open error: %s\n", stafile);
			assert(0);
		}
		fgets(aline, 80, fp_sta);
		while (strncmp(aline + 33, "      ", 6) != 0) {
			nstr = nstr + 1;
			if (nstr > maxlst) {
				printf("  Error:  too many stations in list.. aborting\n");
				fprintf(fp_log,
						"  Error:  too many stations in list.. aborting\n");
				assert(0);
			}

// c
// c----In this version, sty and stx are in meters from the origin - Note that
// c     the y coordinate is in the first column!  slat and slon are lats
// c     and lons in decimal degrees
// c
		}
		printf("  %d stations read in\n", nstr);
		fprintf(fp_log, "  %d stations read in\n", nstr);
		fclose(fp_sta);
	}

	int nstr2 = 2 * nstr;
	int nph = 2;

//c---read corrections
	FILE *fp_fmd = fopen(nmodfil, "rb");
	if (!fp_fmd) {
		printf("file open error: %s\n", nmodfil);
		assert(0);
	}
	int nvar = 0;
	fread(&nvar, sizeof(nvar), 1, fp_fmd);
	if (nvar > maxmbl) {
		printf(
				"  Error:  nvar is greater than maximum: nvar=%12d maxmbl=%12d\n",
				nvar, maxmbl);
		fprintf(fp_log,
				"  Error:  nvar is greater than maximum: nvar=%12d maxmbl=%12d\n",
				nvar, maxmbl);
		assert(0);
	}
	fread(vd, sizeof(vd[0]), nvar, fp_fmd);
	fread(jndx, sizeof(jndx[0]), nvar, fp_fmd);
	for (int i = 0; i < nvar; i++) {
		jndx[i]--;
	}
	fclose(fp_fmd);

	printf("  %10d perturbations read in\n", nvar);
	fprintf(fp_log, "  %10d perturbations read in\n", nvar);

//c---scale the corrections if required
	if (ipscflg == 1) {
		for (int i = 0; i < nvar; i++) {
			vd[i] *= pertscl;
		}
	}

	int isetp = 0, isets = 0;
	int maxvarc = nxc * nyc * nzc;	
	int maxvarc2 = maxvarc * 2;
	if (ido1d == 1) {
		maxvarc = nzc - 1;
		maxvarc2 = maxvarc * 2;
	}

//c---expand
	for (int i = 0; i < maxvarc2; i++) {
		du[i] = 0.0;
	}
	float dupmax = FLT_MIN, dupmin = FLT_MAX;
	float dusmax = FLT_MIN, dusmin = FLT_MAX;
	int istcor = 0;
	for (int i = 0; i < nvar; i++) {
		if (jndx[i] < maxvarc2) {
			if (jndx[i] < maxvarc) {
				if (isetp == 0) {
					dupmax = vd[i];
					dupmin = vd[i];
					isetp = 1;
				} else {
					if (vd[i] > dupmax)
						dupmax = vd[i];
					if (vd[i] < dupmin)
						dupmin = vd[i];
				}
			} else {
				if (isets == 0) {
					dusmax = vd[i];
					dusmin = vd[i];
					isets = 1;
				} else {
					if (vd[i] > dusmax)
						dusmax = vd[i];
					if (vd[i] < dusmin)
						dusmin = vd[i];
				}
			}
			du[jndx[i]] = vd[i];
// c---update station corrections
		} else if (istcor == 1 && (jndx[i] < maxvarc2 + nstr2)) {
			int nst = jndx[i] - nxyzc2;
			if (nst > nstr) {
				nst -= nstr;
				scor[nst] += vd[i];
			} else {
				pcor[nst] += vd[i];
			}
		}
	}
	printf("  Original perturbation extrema from runlsqr file:\n");
	if (isetp == 1)
		printf("  dupmax, dupmin =   %12.8E %12.8E\n", dupmax, dupmin);
	if (isets == 1)
		printf("  dusmax, dusmin =   %12.8E %12.8E\n", dusmax, dusmin);
	fprintf(fp_log, " Original perturbations from runlsqr file:\n");
	if (isetp == 1)
		fprintf(fp_log, "  dupmax, dupmin =   %12.8E %12.8E\n", dupmax, dupmin);
	if (isets == 1)
		fprintf(fp_log, "  dusmax, dusmin =   %12.8E %12.8E\n", dusmax, dusmin);

//c---update station corrections
	if (istacor == 1) {
		FILE *fp_nst = fopen(nstafil, "w");
		if (!fp_nst) {
			printf("file open error: %s\n", nstafil);
			assert(0);
		}
		for (int i = 0; i < nstr; i++) {
			strcpy(aline, bline[i]);
			fprintf(fp_nst, "%62s %f %f\n", aline, pcor[i], scor[i]);
		}
		fprintf(fp_nst, "\n");
		fclose(fp_nst);
	}

//c---read in old velocity model
	char hdr[120];
	FILE *fp_cor = fopen(oldvfil, "r");
	if (!fp_cor) {
		printf("file open error: %s\n", nstafil);
		assert(0);
	}
	struct vhead headin;
	fread(&headin, sizeof(char), nhbyte, fp_cor);
	
	char head[5], type[5], syst[5];
	char quant[5];
	char flatten[5];
	char hcomm[100];
	{
		char *offset = headin.header;
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
		sscanf(offset, "%99s", hcomm);
	}
	fxs = headin.fxs;
	fys = headin.fys;
	fzs = headin.fzs;
	clath = headin.clat;
	clonh = headin.clon;
	czh = headin.cz;
	axo = headin.x0;
	ayo = headin.y0;
	azo = headin.z0;
	dx = headin.dx;
	dy = headin.dy;
	dz = headin.dz;
	azh = headin.az;
	nxh = headin.nx;
	nyh = headin.ny;
	nzh = headin.nz;

// c---see if this is a valid header
	if (strcmp(head, "HEAD") != 0) {
		printf(
				"File does not contain valid header...attempting headerless read \n");
		fprintf(fp_log,
				"File does not contain valid header...attempting headerless read \n");
		// len_rec = 4 * nxyzc2;
		rewind(fp_cor);
		if (!fp_cor) {
			printf("(Error in c2f.c)read file error.\n");
			assert(fp_cor);
		}
		fread(vo, sizeof(vo[0]), nxyzc2, fp_cor);

//c---set trial header values
		strcpy(head, "HEAD");
		strcpy(syst, "CART");
		strcpy(quant, "BMOD");
		strcpy(flatten, "NOFL");
		get_vars(fp_spc, "clat ", pval, &len, &ierr);
		if (ierr == 0)
			sscanf(pval, "%lf", &clath);
		get_vars(fp_spc, "clon ", pval, &len, &ierr);
		if (ierr == 0)
			sscanf(pval, "%lf", &clonh);
		get_vars(fp_spc, "cz ", pval, &len, &ierr);
		if (ierr == 0)
			sscanf(pval, "%lf", &czh);
		get_vars(fp_spc, "azmod ", pval, &len, &ierr);
		if (ierr == 0)
			sscanf(pval, "%f", &azh);
		get_vars(fp_spc, "h ", pval, &len, &ierr);
		if (ierr == 0)
			sscanf(pval, "%lf", &h);
		get_vars(fp_spc, "x0 ", pval, &len, &ierr);
		if (ierr == 0)
			sscanf(pval, "%lf", &axo);
		get_vars(fp_spc, "y0 ", pval, &len, &ierr);
		if (ierr == 0)
			sscanf(pval, "%lf", &ayo);
		get_vars(fp_spc, "z0 ", pval, &len, &ierr);
		if (ierr == 0)
			sscanf(pval, "%lf", &azo);
		dx = h;
		dy = h;
		dz = h;
		nxh = nxc;
		nyh = nyc;
		nzh = nzc;

// c-----Grid specs
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
	} else {
		if (strcmp(type, "CORS") != 0) {
			printf("  WARNING: input mesh does not appear to be COARSE: %s\n",
					type);
			fprintf(fp_log,
					"  WARNING: input mesh does not appear to be COARSE: %s\n",
					type);
		}
		if (strcmp(quant, "BMOD") != 0) {
			printf("  WARNING: file does not appear to be a valid type: %s\n",
					quant);
			fprintf(fp_log,
					"  WARNING: file does not appear to be a valid type: %s\n",
					quant);
		}

		printf("  Reading in Coarse Mesh ...");
		// int len_grd = 4 * (nxc + nyc + nzc - 3);
		// len_rec = len_head + len_grd + 4 * nxyzc2;

		fread(igridx, sizeof(igridx[0]), nxc - 1, fp_cor);
		fread(igridy, sizeof(igridy[0]), nyc - 1, fp_cor);
		fread(igridz, sizeof(igridz[0]), nzc - 1, fp_cor);
		fread(vo, sizeof(vo[0]), nxyzc2, fp_cor);
	}
	fclose(fp_cor);
	printf(".. Done \n");

	for (int i = 0; i < nxyzc2; i++) {
		if (vo[i] < 0) {
			printf("  Error:  %d %f\n", i, vo[i]);
		}
	}

// c----run through the corrections and make sure that none exceed the cutoff
// c	Vo*(1-dvperc) < Vn < Vo*(1+dvperc).
// c    Note that the perturbations are slowness du = 1/Vn - 1/Vo so
// c	1/Vn = du + 1/Vo  or Vn = 1./(du + 1/Vo).  Also, note that
// c	du = -dv/Vn.
// c

// untested
	if (limitu == 1) {
		float dvm1 = 1.f - dvperc;
		float dvm2 = 1.f + dvperc;
		for (int i = 0; i < maxvarc; i++) {
			int ipv = i;
			if (ido1d == 1)
				ipv = nxyc * i;
			float voi = vo[ipv];
			float vnew = 1.f / (du[i] + 1.f / voi);
			float dvn = vnew - voi;
			float vnewmin = voi * dvm1;
			float vnewmax = voi * dvm2;
			if (vnew > vnewmax) {
				dvn = vnewmax - voi;
				du[i] = 1.f / vnewmax - 1.f / voi;
			}
			if (vnew < vnewmin) {
				dvn = vnewmin - voi;
				du[i] = 1.f / vnewmin - 1.f / voi;
			}
			vn[i] = voi + dvn;
		}
// c	Vs or Vp/Vs section
// c	  do i = nxyzc+1, nxyzc2

		for (int i = maxvarc; i < maxvarc2; i++) {
			int isv = i;
			if (ido1d == 1)
				isv = nxyzc + nxyc * (i - maxvarc);
			float voi = vo[isv];
			float vnewmin = voi * dvm1;
			float vnewmax = voi * dvm2;
// c    Case of du = d(Vp/Vs)
			if (ivpvs == 1) {
				float rold = vo[isv - nxyzc] / voi;
				float rnew = du[i] + rold;
				float vnew = vn[i - maxvarc] / rnew;
				if (vnew > vnewmax) {
					rnew = vn[i - maxvarc] / vnewmax;
					du[i] = rnew - rold;
				}
				if (vnew < vnewmin) {
					rnew = vn[i - maxvarc] / vnewmin;
					du[i] = rnew - rold;
				}
			} else {
//c     Case of du = dVS
				float vnew = 1.f / (du[i] + 1.f / voi);
				if (vnew > vnewmax) {
					du[i] = 1.f / vnewmax - 1.f / voi;
				}
				if (vnew < vnewmin) {
					du[i] = 1.f / vnewmin - 1.f / voi;
				}
			}
		}
	}

// c---modification 5/1/03 SWR: we presume that the changes outside the
// c	model box are 0, so to reflect this we average in each case by the total
// c	volume of the averaging box and not just that part that is within
// c	the box.  The effect of using only parts within the box is that
// c	edge effects are not smoothed as much as the remainder.   This is an
// c	important consideration at the base of a teleseismic inversion.
// c
// c
// c	For now, skip smoothing for 1D
	int ioff = 0;

	if (ido1d == 1) {
		for (int i = 0; i < maxvarc2; i++) {
			ds[i] = du[i];
		}
	} else {
		int kd = 2 * kma + 1;
		int jd = 2 * jma + 1;
		int id = 2 * ima + 1;
		int ntot = kd * jd * id;
// c  Smooth the perturbations
		for (int iii = 0; iii < nsmooth; iii++) {
			for (int n = 0; n < nph; n++) {
				ioff = nxyzc * n;
				for (int k = 0; k < nzc; k++) {
					int nk = nxyc * k + ioff;
					int kmin = k - kma;
					int kmax = k + kma;
					if (kmin < 0)
						kmin = 0;
					if (kmax > nzc - 1)
						kmax = nzc - 1;
					// int numz = kmax - kmin;
					for (int j = 0; j < nyc; j++) {
						int nkj = nk + nxc * j;
						int jmin = j - jma;
						int jmax = j + jma;
						if (jmin < 0)
							jmin = 0;
						if (jmax > nyc - 1)
							jmax = nyc - 1;
						// int numyz = (jmax - jmin) * numz;
						int i = 0;
						int iiiii = nkj + i;
						ds[iiiii] = 0;
						int imin = 0;
						int imax = i + ima;
						if (imax > nxc - 1)
							imax = nxc - 1;
						// int num = (imax - imin) * numyz;
						for (int kk = kmin; kk <= kmax; kk++) {
							int nnkk = nxc * nyc * kk + ioff;
							for (int jj = jmin; jj <= jmax; jj++) {
								int nnkkjj = nnkk + nxc * jj;
								//printf("index=%d, jj=%d, kk=%d\n", nnkkjj, jj, kk);
								for (int ii = imin; ii <= imax; ii++) {
									ds[iiiii] += du[nnkkjj + ii];
								}							
							}
							//printf("\n");
						}
						float vlast = ds[iiiii];
						ds[iiiii] /= ntot;
// c              COMPUTE FOR REMAINING i VALUES
						for (int i = 1; i < nxc; i++) {
							int imin = i - ima;
							int imax = i + ima;
							if (imin < 0)
								imin = 0;
							if (imax > nxc - 1)
								imax = nxc - 1;
							float walll = 0;
							float wallr = 0;
							// num = (imax - imin) * numyz;
							if (imin > 0) {
								for (int kk = kmin; kk <= kmax; kk++) {
									int nnkk = nxc * nyc * kk + ioff;
									for (int jj = jmin; jj <= jmax; jj++) {
										walll += du[nnkk + nxc * jj + imin - 1];
									}
								}
							}
							if (i + ima < nxc) {
								for (int kk = kmin; kk <= kmax; kk++) {
									int nnkk = nxc * nyc * kk + ioff;
									for (int jj = jmin; jj <= jmax; jj++) {
										wallr += du[nnkk + nxc * jj + imax];
									}
								}
							}
							iiiii++;
							ds[iiiii] = vlast - walll + wallr;
							vlast = ds[iiiii];
							ds[iiiii] /= ntot;
						}
					}
				}

				float dumin = FLT_MAX;
				float dumax = FLT_MIN;
				for (int k = 0; k < nxyzc; k++) {
					int i = k + ioff;
					du[i] = ds[i];
					if (k == 0) {
						dumin = du[i];
						dumax = du[i];
					} else {
						if (du[i] < dumin)
							dumin = du[i];
						if (du[i] > dumax)
							dumax = du[i];
					}
				}
				printf("  nsmooth, n, dumax, dumin = %12d%12d %15.8E %15.8E\n",
						iii + 1, n + 1, dumax, dumin);
				fprintf(fp_log,
						" nsmooth, n, dumax, dumin = %12d%12d %15.8E %15.8E\n",
						iii + 1, n + 1, dumax, dumin);
			}
		}
	}

// c----redefine vo to be Vp/Vs if necessary (perturbation is d(vp/vs) if ivpvs is 1).
	if (ivpvs == 1) {
		for (int i = 0; i < nxyzc; i++) {
			vo[i + nxyzc] = vo[i] / vo[i + nxyzc];
		}
	}

// c-----add corrections to existing model
// c     TEST FOR 2D MODELS
// c--NB:  this part of the code (2D) has NOT BEEN TESTED!
	int ix2d = 0, iy2d = 0; // , iz2d = 0;
	if (nxc == 1) {
		ix2d = 1;
		printf("  2d model encountered (nx=1)\n");
		printf("  WARNING: THIS PART OF THE CODE IS UNTESTED!\n");
		for (int nn = 0; nn < 2; nn++) {
			for (int k = nzc - 1; k >= 0; k--) {
				int nk = nxc * nyc * k;
				for (int j = nyc - 1; j >= 0; j--) {
					int iiiii = nk + nxc * j + 1;
					ds[iiiii] = ds[nyc * k + j + 1];
					iiiii++;
					ds[iiiii] = ds[nyc * k + j + 1];
				}
			}
		}
	}
	if (nyc == 1) {
		iy2d = 1;
		printf("  2d model encountered (ny=1)\n");
		printf("  WARNING: THIS PART OF THE CODE IS UNTESTED!\n");
		for (int nn = 0; nn < 2; nn++) {
			for (int k = nzc - 1; k >= 0; k--) {
				int nk = nxc * nyc * k;
				for (int i = nyc - 1; i >= 0; i--) {
					int iiiii = nk + i;
					ds[iiiii] = ds[nxc * k + i];
					iiiii += nxc;
					ds[iiiii] = ds[nxc * k + i];
				}
			}
		}
	}
	if (nzc == 1) {
		iy2d = 1;
		printf("  2d model encountered (nz=1)\n");
		printf("  WARNING: THIS PART OF THE CODE IS UNTESTED!\n");
		for (int nn = 0; nn < 2; nn++) {
			for (int j = 0; j < nyc; j++) {
				for (int i = 0; i < nxc; i++) {
					int iiiii = nxc * nyc + nxc * j + i;
					ds[iiiii] = ds[nx * j + i];
				}
			}
		}
	}
// c--end of 2D test

// c---add Vp perturbation to model
	float dv_min = FLT_MAX, dv_max = FLT_MIN;
	float v_min_new = FLT_MAX, v_max_new = FLT_MIN;
	float v_min_old = FLT_MAX, v_max_old = FLT_MIN;
	float ddu = 0, vold = 0, vnew = 0, dv = 0;
	for (int i = 0; i < nxyzc; i++) {
		int ivar = i;
		if (ido1d == 1)
			ivar = i / nxyc;
		ddu = ds[ivar];
		vold = vo[i];
		vnew = 1.f / (ddu + (1.f / vold));
		dv = vnew - vold;
		if (i == 0) {
			dv_min = dv;
			dv_max = dv;
			v_min_new = vnew;
			v_max_new = vnew;
			v_min_old = vold;
			v_max_old = vold;
		} else {
			if (dv > dv_max)
				dv_max = dv;
			if (dv < dv_min)
				dv_min = dv;
			if (vnew > v_max_new)
				v_max_new = vnew;
			if (vnew < v_min_new)
				v_min_new = vnew;
			if (vold > v_max_old)
				v_max_old = vold;
			if (vold < v_min_old)
				v_min_old = vold;
		}
		vn[i] = vnew;
	}
	printf("  dvpmax, dvpmin = %15.8E %15.8E\n", dv_max, dv_min);
	printf("  Old Vp max, New Vp max = %12.8f %15.8f\n", v_max_old, v_max_new);
	printf("  Old Vp min, New Vp min = %12.8f %15.8f\n", v_min_old, v_min_new);
	fprintf(fp_log, "  dvpmax, dvpmin = %15.8E %15.8E\n", dv_max, dv_min);
	fprintf(fp_log, "  Old Vp max, New Vp max = %12.8f %15.8f\n", v_max_old,
			v_max_new);
	fprintf(fp_log, "  Old Vp min, New Vp min = %12.8f %15.8f\n", v_min_old,
			v_min_new);

// c---add Vs perturbation to model - recall that vo was changed to Vp/Vs above if ivpvs = 1
	dv = 0, ddu = 0, vold = 0;
	dv_min = FLT_MAX, dv_max = FLT_MIN;
	v_min_new = FLT_MAX, v_max_new = FLT_MIN;
	v_min_old = FLT_MAX, v_max_old = FLT_MIN;
	ioff = nxyzc;
	for (int k = 0; k < nxyzc; k++) {
		int i = k + ioff;
		int ivar = i;
		if (ido1d == 1)
			ivar = k / nxyc + nzc;
		ddu = ds[ivar];
		vold = vo[i];
		if (ivpvs == 1) {
			vnew = vold + ddu;
			dv = ddu;
		} else {
			vnew = 1 / (ddu + (1 / vold));
			dv = vnew - vold;
		}

		if (k == 0) {
			dv_min = dv;
			dv_max = dv;
			v_min_new = vnew;
			v_max_new = vnew;
			v_min_old = vold;
			v_max_old = vold;
		} else {
			if (dv > dv_max)
				dv_max = dv;
			if (dv < dv_min)
				dv_min = dv;
			if (vnew > v_max_new)
				v_max_new = vnew;
			if (vnew < v_min_new)
				v_min_new = vnew;
			if (vold > v_max_old)
				v_max_old = vold;
			if (vold < v_min_old)
				v_min_old = vold;
		}
		vn[i] = vnew;
	}

	if (ivpvs == 1) {
		printf("  d(vp/vs)max, d(vp/vs)smin = %15.8E %15.8E\n", dv_max, dv_min);
		printf("  Old Vp/Vs max, New Vp/Vs max = %12.8f %15.8f\n", v_max_old,
				v_max_new);
		printf("  Old Vp/Vs min, New Vp/Vs min = %12.8f %15.8f\n", v_min_old,
				v_min_new);
		fprintf(fp_log, "  d(vp/vs)max, d(vp/vs)smin = %15.8E %15.8E\n", dv_max,
				dv_min);
		fprintf(fp_log, "  Old Vp/Vs max, New Vp/Vs max = %12.8f %15.8f\n",
				v_max_old, v_max_new);
		fprintf(fp_log, "  Old Vp/Vs min, New Vp/Vs min = %12.8f %15.8f\n",
				v_min_old, v_min_new);
	} else {
		printf("  dvsmax, dvsmin = %15.8E %15.8E\n", dv_max, dv_min);
		printf("  Old Vs max, New Vs max = %12.8f %15.8f\n", v_max_old,
				v_max_new);
		printf("  Old Vs min, New Vs min = %12.8f %15.8f\n", v_min_old,
				v_min_new);
		fprintf(fp_log, "  dvsmax, dvsmin = %15.8E %15.8E\n", dv_max, dv_min);
		fprintf(fp_log, "  Old Vs max, New Vs max = %12.78 %15.8f\n", v_max_old,
				v_max_new);
		fprintf(fp_log, "  Old Vs min, New Vs min = %12.8f %15.8f\n", v_min_old,
				v_min_new);
	}
	fclose(fp_log);

// c---convert Vp/Vs back to Vs if necessary
	if (ivpvs==1) {
	  for(int k=0;k<nxyzc;k++) {
	    int i = k + nxyzc;
	    vn[i] = vn[k]/vn[i];
	  }
	}

// c     WRITE NEW MODEL (CHECK TO HANDLE 2D FIRST)
	if (ix2d==1) {
		int nxs = 0;
		for(int k=0; k<nzc;k++) {
			int nk = nxs*nyc*k;
			for(int j=0;j<nyc;j++) {
				int iiiii = nk+nxs*j;
				vn[iiiii] = vn[2*nyc*k+2*j];
			}
		}
	}
	if (iy2d==1) {
		int nys = 0;
		for(int k=0; k<nzc;k++) {
			int nk = nxc*nys*k;
			for(int i=0;i<nxc;i++){
				int iiiii = nk+i;
				vn[iiiii] = vn[nxc*2*k+i];
			}
		}
	}
	// if (iz2d==1) nzs = 1;

	trim(fmodfil);
	char vpfile[strlen(fmodfil)+6], vsfile[strlen(fmodfil)+6];
	sprintf(vpfile,"%s.pvel", fmodfil);
	sprintf(vsfile,"%s.svel", fmodfil);
	struct vhead headout;
	headout.fxs = fxs;
	headout.fys = fys;
	headout.fzs = fzs;
	headout.clat = clath;
	headout.clon = clonh;
	headout.cz = czh;
	headout.x0 = axo;
	headout.y0 = ayo;
	headout.z0 = azo;
	headout.dx = dx;
	headout.dy = dy;
	headout.dz = dz;
	headout.az = azh;
	headout.nx = nxh;
	headout.ny = nyh;
	headout.nz = nzh;

	sprintf(hcomm, "Output from makenewmod.c Version %s using %s", VERSION, oldvfil);
	hdr_appender(hdr, nhbyte,head,type,syst,quant,flatten,hcomm);
	strncpy(headout.header, hdr, 120);
	FILE *fp_nmd=fopen(fmodfil, "wb");
	fwrite(&headout, sizeof(char), nhbyte, fp_nmd);
	fwrite(igridx, sizeof(igridx[0]), nxc - 1, fp_nmd);
	fwrite(igridy, sizeof(igridy[0]), nyc - 1, fp_nmd);
	fwrite(igridz, sizeof(igridz[0]), nzc - 1, fp_nmd);
	fwrite(vn, sizeof(vn[0]), nxyzc2, fp_nmd);
	fclose(fp_nmd);

	hdr_appender(hdr, nhbyte,head,type,syst,"VPMD",flatten,hcomm);
	strncpy(headout.header, hdr, 120);
	FILE *fp_vp=fopen(vpfile, "wb");
	fwrite(&headout, sizeof(char), nhbyte, fp_vp);
	fwrite(igridx, sizeof(igridx[0]), nxc - 1, fp_vp);
	fwrite(igridy, sizeof(igridy[0]), nyc - 1, fp_vp);
	fwrite(igridz, sizeof(igridz[0]), nzc - 1, fp_vp);
	fwrite(vn, sizeof(vn[0]), nxyzc, fp_vp);
	fclose(fp_vp);

	hdr_appender(hdr, nhbyte,head,type,syst,"VSMD",flatten,hcomm);
	strncpy(headout.header, hdr, 120);
	FILE *fp_vs=fopen(vsfile, "wb");
	fwrite(&headout, sizeof(char), nhbyte, fp_vs);
	fwrite(igridx, sizeof(igridx[0]), nxc - 1, fp_vs);
	fwrite(igridy, sizeof(igridy[0]), nyc - 1, fp_vs);
	fwrite(igridz, sizeof(igridz[0]), nzc - 1, fp_vs);
	for(int i=nxyzc;i<nxyzc2;i++) {
		fwrite(vn, sizeof(vn[0]), 1, fp_vs);
	}
	fclose(fp_vs);

// c---compute a 1D average for the elements that were hit
	FILE *fp_1dm=fopen("new1d.mod", "w");
	get_vars(fp_spc, "z0 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &z0);
	get_vars(fp_spc, "h ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &h);
	fclose(fp_spc);

	float gz[nzcm];
	gz[0]=z0;
	for(int i=1;i<nzc;i++) {
		gz[i]=gz[i-1]+h*igridz[i-1];
	}

	for(int i=1;i<nzc*2;i++) {
		du[i]=0;
		igridz[i]=0;
	}

	for(int i=0;i<nvar;i++) {
		int jnd=jndx[i];
		if(jnd<nxyzc2) {
			int iph=jnd/nxyzc;
			int kz=(jnd-iph*nxyzc)/nxyc;
			int iz=kz+iph*nzc;
			du[iz]+=vn[jnd];
			igridz[iz]++;
		}
	}

	for(int i=0;i<nzc;i++) {
		if(igridz[i]>0) 
			du[i]/=igridz[i];
		if(igridz[i+nzc]>0)
			du[i+nzc]/=igridz[i+nzc];
		if(du[i+nzc]>-0.001 && du[i+nzc]<0.001)
			du[i+nzc]=du[i]/1.78f;
		fprintf(fp_1dm, "%7.2f%10.5f%10.5f I\n", gz[i], du[i], du[i+nzc]);
	}
	fclose(fp_1dm);
}
