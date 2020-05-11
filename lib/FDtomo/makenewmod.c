// c
// c	Program to add perturbations to a wavespeed model to that model
// c	and output a new model.  The perturbations can be clipped and smoothed
// c	via moving average windows if desired.
// c
// c	Existing and output models are a Coarse Mesh grid wavespeed file with header.
// c	Pertubation file is generated by runlsqr.f
// c
// c	VERSION 2004.0909
// c
// c       Input File attachments:
// c
// c       fn      Name   Uname   U#  Comments
// c	36     nmodfil lunfmd  45  Perturbation file (from runf04 or runlsqr)
// c	 5     oldvfil luncor  12  Old coarse scale model
// c	 1     stafile lunsta   8  Station list (if istacor = 1)
// c
// c       Output File attachments:
// c
// c       fn      Name   Uname   U#  Comments
// c                      lout     2  makenewmod.log Log File
// c	38     fmodfil lunnmd  47  P&S new model
// c	39     modpfil lunnmp  48  P new model (only filenum used)
// c	40     modsfil lunnms  49  S new model (only filenum used)
// c	       nstafil lunnst  15  New station list (if istacor = 1)
// c
// c	Note on limits and smoothing:
// c
// c	The user can specify a maximum percentage change in wavespeed for the
// c	model by setting the parameter dvperc = |Vn - Vo|/Vo, where Vo
// c	is the existing wavespeed.  In general, this means that 
// c	the new wavespeed (Vn) is limited relative to the existing wavespeed Vo by
// c	Vo*(1-dvperc) < Vn < Vo*(1+dvperc).  So, for each perturbation we check
// c	that Vn is within these limits, and if not we reset it to those limits.
// c
// c	Required variables:
// c
// c	nxc, nyc, nzc	Number of grid points in x, y, z directions
// c
// c	Optionally  specified variables:
// c
// c	pscflg
// c		Flag to allow scaling of perturbations at the outset.  
// c		=1 to scale; =0 not to scale
// c		Default = 0
// c
// c	pertscl
// c		fraction to scale the perturbations at the outset (if pscflg = 1)
// c	    	applied as vd(i) = vd(i)*pertscl
// c
// c	dvperc 
// c		Maximum percentage change allowd for percentage of dv
// c	 	Default = 0.05
// c
// c	limitu
// c		Flag to apply maximum limits on pertubarions:
// c		= to apply limits; = 0 do not apply.
// c		Default = 0
// c
// c	mavx, mavy, mavz 
// c		Moving average window length in the (x, y, z) direction. 
// c		Default = 3, 3, 3
// c
// c	nsmooth
// c		Number of iterations over which the moving average window is applied.
// c		Default = 2
// c
// c	ivpvs   
// c		if = 1, then the pertubations to the S part of the model are
// c		d(Vp/Vs) rather than dVs
// c		Default = 1
// c
// c	istacor 
// c		if = 1, then read in and apply station corrections
// c		Default = 0
// c
// c	do1d
// c		if = 1, then add the same perturbation to all wavespeeds at the
// c		same depth (perturbations will be ordered by depth).
// c

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
#include "common/shared_variables.h"
#include "common/read_spec.h"
#include "FDtomo/makenewmod.h"


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
#define VERSION "2018.0625"

// c---du holds all the perturbations.  
// c   ds holds the smoothed perturbations
// c   vo holds the original model
// c   vn holds the new model
float du[nxyzcm2], ds[nxyzcm2], *vn;

// c---vd are the pertubations from runlsqr
float *vd;
int *jndx;

MAKENEWMOD_DATA *makenewmod(Coordinate3D coordinate, velocityModel3D vp_model, velocityModel3D vs_model, RUNLSQR_DATA *RUNLSQR,
							 int numberOfStations, MakenewmodEnv makenewmod_env, CommonEnv common_env) {
	MAKENEWMOD_DATA *MAKENEWMOD = (MAKENEWMOD_DATA *)malloc(sizeof(MAKENEWMOD_DATA));
	//initialize variable
	int nxc = coordinate.mesh.numberOfNode.x;
	int nyc = coordinate.mesh.numberOfNode.y;
	int nzc = coordinate.mesh.numberOfNode.z;

	int *igridx = coordinate.mesh.gridx;
	int *igridy = coordinate.mesh.gridy;
	int *igridz = coordinate.mesh.gridz;
	
	char aline[MAXSTRLEN + 1], bline[maxlst][MAXSTRLEN + 1];
	char fmodfil[MAXSTRLEN + 1];

	int istacor = common_env.istacor;
	int ivpvs = common_env.ivpvs;

	int limitu = makenewmod_env.limitu;
	int mavx = makenewmod_env.mavx;
	int mavy = makenewmod_env.mavy;
	int mavz = makenewmod_env.mavz;
	int nsmooth = makenewmod_env.nsmooth;
	int ipscflg = makenewmod_env.ipscflg;
	int ido1d = makenewmod_env.ido1d;
	float dvperc = makenewmod_env.dvperc;
	float pertscl = makenewmod_env.pertscl;

	char nstafil[MAXSTRLEN + 1];

	vn = (float *)malloc(sizeof(float) * nxyzcm2);

	int nxyc = nxc * nyc;
	int nxyzc = nxyc * nzc;
	int nxyzc2 = nxyzc * 2;

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

	int nstr = 0;
	float pcor[maxlst], scor[maxlst];
	memset(pcor, 0, sizeof(pcor));
	memset(scor, 0, sizeof(scor));
	if (istacor == 1) {
		nstr = numberOfStations;
	}

	int nstr2 = 2 * nstr;
	int nph = 2;
	int nvar = 0;

	nvar = RUNLSQR->n;
	vd = RUNLSQR->x;
	jndx = RUNLSQR->jndx;

	for (int i = 0; i < nvar; i++) {
		jndx[i]--;
	}
	printf("  %10d perturbations read in\n", nvar);

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
	//powqkpowkqpok
	float *vp = vp_model.velocity;
	float *vs = vs_model.velocity;
	for (int i = 0; i < nxyzc; i++) {
		if (vp[i] < 0) {
			printf("  Error:  %d %f\n", i, vp[i]);
			assert(0);
		}
		if (vs[i] < 0) {
			printf("  Error:  %d %f\n", i, vs[i]);
			assert(0);
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
			float voi = vp[ipv];
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
		for (int k = 0; k < maxvarc; k++) {
			int ksv = k;
			int i = k + maxvarc;
			int isv = i;
			if (ido1d == 1){
				ksv = nxyc * k;
				isv = nxyzc + nxyc * (i - maxvarc);
			}
			float voi = vs[k];
			float vnewmin = voi * dvm1;
			float vnewmax = voi * dvm2;
// c    Case of du = d(Vp/Vs)
			if (ivpvs == 1) {
				float rold = vp[ksv] / voi;
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
			}
		}
	}

// c----redefine vo to be Vp/Vs if necessary (perturbation is d(vp/vs) if ivpvs is 1).
	if (ivpvs == 1) {
		for (int i = 0; i < nxyzc; i++) {
			vs[i] = vp[i] / vs[i];
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
					//ds[iiiii] = ds[nx * j + i];
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
		vold = vp[i];
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

// c---add Vs perturbation to model - recall that vo was changed to Vp/Vs above if ivpvs = 1
	dv = 0, ddu = 0, vold = 0;
	dv_min = FLT_MAX, dv_max = FLT_MIN;
	v_min_new = FLT_MAX, v_max_new = FLT_MIN;
	v_min_old = FLT_MAX, v_max_old = FLT_MIN;
	ioff = nxyzc;
	for (int k = 0; k < nxyzc; k++) {
		int s_index = k;
		int i = k + ioff;
		int ivar = i;
		if (ido1d == 1){
			ivar = k / nxyc + nzc;
			s_index = k / nxyc;
		}
		ddu = ds[ivar];
		vold = vs[s_index];
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
	} else {
		printf("  dvsmax, dvsmin = %15.8E %15.8E\n", dv_max, dv_min);
		printf("  Old Vs max, New Vs max = %12.8f %15.8f\n", v_max_old,
				v_max_new);
		printf("  Old Vs min, New Vs min = %12.8f %15.8f\n", v_min_old,
				v_min_new);
	}

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
	memcpy(MAKENEWMOD->igridx, igridx, sizeof(igridx[0]) * (nxc - 1));
	memcpy(MAKENEWMOD->igridy, igridy, sizeof(igridy[0]) * (nyc - 1));
	memcpy(MAKENEWMOD->igridz, igridz, sizeof(igridz[0]) * (nzc - 1));
	MAKENEWMOD->vn = vn;

// c---compute a 1D average for the elements that were hit
	FILE *fp_1dm=fopen("new1d.mod", "w");
	double z0 = coordinate.origin.z;
	double h = coordinate.space.z;

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
	
	return MAKENEWMOD;
}

int OUTPUT_MAKENEWMOD(MAKENEWMOD_DATA *MAKENEWMOD, SPEC spec){
	FILE *fp_nmd;
	fp_nmd = fopen(spec.fmodfil, "wb");
	if (!fp_nmd) {
		printf("file create error: %s\n", spec.fmodfil);
		assert(0);
	}

	int nxc = spec.grid.nxc, nyc = spec.grid.nyc, nzc = spec.grid.nzc;
	int nxyc = nxc * nyc;
	int nxyzc = nxyc * nzc;
	int nxyzc2 = nxyzc * 2;
	fwrite(&MAKENEWMOD->head, sizeof(char), nhbyte, fp_nmd);
	fwrite(MAKENEWMOD->igridx, sizeof(MAKENEWMOD->igridx[0]), nxc - 1, fp_nmd);
	fwrite(MAKENEWMOD->igridy, sizeof(MAKENEWMOD->igridy[0]), nyc - 1, fp_nmd);
	fwrite(MAKENEWMOD->igridz, sizeof(MAKENEWMOD->igridz[0]), nzc - 1, fp_nmd);
	fwrite(MAKENEWMOD->vn, sizeof(MAKENEWMOD->vn[0]), nxyzc2, fp_nmd);

	return 0;
}

int LOG_MAKENEWMOD(SPEC spec){
	char logfile[80 + 1];
	sprintf(logfile, "sphfdloc.log%d", spec.ittnum);
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
			spec.spec_file);
	fprintf(fp_log, "  \n");
	{
		char tmp[MAXSTRLEN];
		dtoa(tmp, spec.ittnum, 18);
		fprintf(fp_log, " Iteration counter:          %s     \n", tmp);
	}
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of X fine grid nodes: %12d\n", spec.grid.nxc);
	fprintf(fp_log, "  Number of Y fine grid nodes: %12d\n", spec.grid.nyc);
	fprintf(fp_log, "  Number of Z fine grid nodes: %12d\n", spec.grid.nzc);
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of X fine grid nodes: %12d\n", spec.grid.nxc);
	fprintf(fp_log, "  Number of Y fine grid nodes: %12d\n", spec.grid.nyc);
	fprintf(fp_log, "  Number of Z fine grid nodes: %12d\n", spec.grid.nzc);
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Moving Window Length in X : %12d\n", spec.mavx);
	fprintf(fp_log, "  Moving Window Length in Y : %12d\n", spec.mavy);
	fprintf(fp_log, "  Moving Window Length in Z : %12d\n", spec.mavz);
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of smoothing passes : %12d\n", spec.nsmooth);
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  ivpvs: %12d\n", spec.ivpvs);
	fprintf(fp_log, " \n");

	if (spec.ipscflg == 1) {
		fprintf(fp_log, " Perturbations scaled at the outset by : %f", spec.pertscl);
	} else {
		fprintf(fp_log, " Perturbations not scaled at the outset (pscflg = 0)");
	}
	if (spec.limitu == 1) {
		fprintf(fp_log, " Maximum Wavespeed Percentage change allowed : %f",
				spec.dvperc);
	} else {
		fprintf(fp_log, " No perturbation limits applied (limitu = 0) ");
	}
	if (spec.ido1d == 1) {
		fprintf(fp_log, " Corrections assumed for 1D model (do1d = 1) ");
	} else {
		fprintf(fp_log, " Corrections assumed for 3D model (do1d = 0) ");
	}

	fprintf(fp_log, " Input file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Arr. Time Wavespeed model: %s\n", spec.oldvfil);
	fprintf(fp_log, " Perturbation file: %s\n", spec.nmodfil);
	fprintf(fp_log, " Existing Station List: %s\n", spec.stafile);
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Output file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " New Wavespeed Model: %s\n", spec.fmodfil);
	fprintf(fp_log, " New Station List : %s\n", spec.nstafil);
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Total Number of coarse grid nodes:%13d\n", spec.grid.nxc * spec.grid.nyc * spec.grid.nzc);
	fprintf(fp_log, " \n");
	fprintf(fp_log,
			" *************************************************************** \n");

	fclose(fp_log);
	return 0;
}

MakenewmodEnv setMakeNewmodEnv(char *spec_file){
	MakenewmodEnv makenewmod_env;
	//default settings
	makenewmod_env.mavx = 3;
	makenewmod_env.mavy = 3;
	makenewmod_env.mavz = 3;
	makenewmod_env.nsmooth = 2;
	makenewmod_env.limitu = 0;
	makenewmod_env.ipscflg = 0;
	makenewmod_env.dvperc = 0.05;
	makenewmod_env.pertscl = 1.0;
	makenewmod_env.ido1d = 0;
	return makenewmod_env;
}

void setMakeNewmodVariables(MakenewmodEnv *makenewmod_env, char *spec_file){
	FILE *fp_spc;
    fp_spc = fopen(spec_file, "r");
	if (!fp_spc) {
		printf("(Error in read_spec.c)read fp_spc file error.\n");
		assert(0);
	}
    
    int len, ierr;
    char pval[MAXSTRLEN + 1];

	get_vars(fp_spc, "mavx ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &makenewmod_env->mavx);
	get_vars(fp_spc, "mavy ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &makenewmod_env->mavy);
	get_vars(fp_spc, "mavz ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &makenewmod_env->mavz);
	get_vars(fp_spc, "nsmooth ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &makenewmod_env->nsmooth);
	get_vars(fp_spc, "limitu ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &makenewmod_env->limitu);
	get_vars(fp_spc, "ipscflg ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &makenewmod_env->ipscflg);
	get_vars(fp_spc, "ipscflg ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &makenewmod_env->ipscflg);
	get_vars(fp_spc, "pertscl ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &makenewmod_env->pertscl);
	get_vars(fp_spc, "do1d ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &makenewmod_env->ido1d);
}
