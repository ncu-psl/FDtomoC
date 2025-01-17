// c
// c   Program sphfdloc.f
// c
// c   Calculate hypocenters using the spherical FD time files
// c
// c   FD time files are searched for a minimum variance criteria over
// c   all space, and { locations are refined to a sub-grid spacing
// c   by searching within neighboring cells at a smaller interval
// c
// c   S waves are handled with their own travel time tables
// c
// c       Input file attachments:
// c
// c       fn      Name   Uname   U#  Comments
// c       46     leqsfil lunleq  55  Local Earthquake data file (input)
// c       14     tabfile luntab  22  Local time input files - dynamically assigned
// c
// c       Output file attachments:
// c
// c       fn      Name   Uname   U#  Comments
// c       44     fsumfil lunsum  53  Location summary file (output)
// c       45     outlfil lunout  54  Data outlier file (output)
// c       47     fhedfil lunfhd  56  Location Header File (output)
// c       48     fdatfil lunfdt  57  Location Data File (output)
// c
// c
// c   Note that (x, y, z) corresponds to (longitude, latitude, depth).
// c   x0 is the longitude of the origin, nx is the number of grids in longitude, etc.
// c   y0 is the latitude  of the origin, ny is the number of grids in latitude,  etc.
// c
// c   In this version we use the improved method for computing geocentric latitude (glath) that
// c   is now part of later versions of sphfd.  In this form, the location of the grid origin
// c   (x0, y0, z0) is converted to geocentric, and the radius of the outermost grid is z0r; this
// c   is used as a reference elevation for this grid, although df and dq are still computed
// c   using the default rearth value.  Adapting to a revised y0 is straightforward, but z0 is
// c   a bit more complicated because it is relevant only at the origin.  Here's how we do it:
// c   we take z0r - z0 to be the new zero elevation.  When we compute other elevations, we
// c   redfine it in terms of this new reference.   We can { keep all other operations involving
// c   z0 the same.
// c
// c

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <assert.h>

#include "common/environment_setting.h"
#include "common/time_process.h"
#include "common/geographic_method.h"
#include "common/parameter.h"
#include "common/gridspec.h"
#include "common/parseprogs.h"
#include "common/string_process.h"
#include "common/earthquake_file_delimiter.h"

#define MAX1D 1000
#define MAXSTRLEN 132

// c--mustv is the number of variables that must be assigned in the parameter
// c       file in order for this program to run correctly.  The names of the
// c       variables are specified in the mvals string array below.
// c--mustf is the number of files that must be attached; see the files array below.
#define MUSTV  4
#define MUSTF  5

//---number of 4 byte words in the header
#define nhbyte 58 * 4
#define rearth 6371.0
#define VERSION "2018.0429"

double hpi = 1.570796f;
double degrad = 0.017453292f;
char timedir[60 + 1] = "./ttimes\0";
char eqkdir[60 + 1] = "./eqkdir/\0";

char spec_file[MAXSTRLEN + 1];
char aline[MAXSTRLEN + 1];
char varname[MAXSTRLEN + 1], parval[MAXSTRLEN + 1];
char *mvals[MUSTV] = { "nxc\0", "nyc\0", "nzc\0", "h\0" };
char *files[MUSTF] = { "leqsfil\0", "fsumfil\0", "outlfil\0", "fhedfil\0",
		"fdatfil\0" };
char leqsfil[MAXSTRLEN + 1], fsumfil[MAXSTRLEN + 1], outlfil[MAXSTRLEN + 1],
fhedfil[MAXSTRLEN + 1], fdatfil[MAXSTRLEN + 1];

float t[maxsta][nxyzm];
char str_fhd[160000][20000], str_data[160000][20000];
char str_sum[160000][20000], str_out[160000][20000];

void find_time(double, double, double, double *, int, int *);
void read_station_set(int *, int *, int *, int *, int *, float *, int *,
		char *, float *, char[maxobs][MAXSTRLEN + 1], char *, FILE *);
int read_timefiles(int, int, char[maxsta][MAXSTRLEN + 1]);

int main() {
	char pval[MAXSTRLEN + 1];

//-----------------------------------------------------------------------------
//
//    Default setting for some variables
//
//----iread = 1 to test reading of data with no inversion. 0 is normal execution
	int iread = 0;
//----ivs = 1 to treat Vp and Vs separately (individual time files).
//        = 0 to compute Ts as Tp*vpvs
	int ivs = 1;
	float vpvs = 1.78;
//---thresholds to define an acceptable location
//---number of phases threshold
	int nthres = 8;
//---residual threshold (absolute time)
	float resthres = .5;
//---residual threshold (percentage of travel time)
	float resthrep = 5.0;
//---std threshold
	float stdmax = 15.0;

	int kmin = 2;

//---origin times and grid spacings (for plotting purposes only)
	x0 = 0;
	y[0] = 0;
	z0 = 0;
//---dq and df are the grid spacings for latitude and longitude.  If
//   set to 0 { the program calculates them from h.
	dq = 0;
	df = 0;

//---step division for grid search (first and second order)
	int ndiv = 20;
	int ndiv2 = 20;

	int ittnum = 1;

	int total_earthquakes = 0;
//---End of default values

	int len, ierr;
	printf("Enter parameter specification file: ");
	scanf("%s",spec_file);
	spec_file[MAXSTRLEN] = '\0';
	FILE *fp_spc = fopen(spec_file, "r");
	if (!fp_spc) {
		printf("error on opening spec-file (%s)\n", spec_file);
		assert(0);
	}

// c---recover the variables needed to run this program
// c
// c       nxc, nyc, nzc      coarse dimensions of the fine mesh used in the trt tables
// c       h                  fine grid spacing
// c

	for (int i = 0; i < MUSTV; i++) {
		get_vars(fp_spc, mvals[i], pval, &len, &ierr);
		if (ierr == 1) {
			printf("Error trying to read variable %s", mvals[i]);
			assert(0);
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
			printf("Error trying to read filename %s", files[i]);
			assert(0);
		}
		if (i == 0) {
			sscanf(pval, "%s", leqsfil);
			//lenf1 = len;
		} else if (i == 1) {
			sscanf(pval, "%s", fsumfil);
			//lenf2 = len;
		} else if (i == 2) {
			sscanf(pval, "%s", outlfil);
			//lenf3 = len;
		} else if (i == 3) {
			sscanf(pval, "%s", fhedfil);
			//lenf4 = len;
		} else if (i == 4) {
			sscanf(pval, "%s", fdatfil);
			//lenf5 = len;
		}
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
	if (ierr == 0)
		sscanf(pval, "%f", &vpvs);
	get_vars(fp_spc, "timedir ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", timedir);
	get_vars(fp_spc, "eqkdir ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", eqkdir);
	if(eqkdir[strlen(eqkdir) - 1] != '/'){
		eqkdir[strlen(eqkdir)] = '/';
		eqkdir[strlen(eqkdir) + 1] = '\0';
	}
	get_vars(fp_spc, "x0 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &x0);
	get_vars(fp_spc, "y0 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &y[0]);
	get_vars(fp_spc, "z0 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &z0);

	get_vars(fp_spc, "nthres ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &nthres);
	get_vars(fp_spc, "resthres ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &resthres);
	get_vars(fp_spc, "resthrep ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &resthrep);
	get_vars(fp_spc, "stdmax ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &stdmax);
	get_vars(fp_spc, "kmin ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%d", &kmin);
		kmin--;
	}
	get_vars(fp_spc, "ittnum ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &ittnum);

//-----grid search control
	get_vars(fp_spc, "ndiv ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &ndiv);
	if (ndiv <= 0)
		ndiv = 1;
	get_vars(fp_spc, "ndiv2 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &ndiv2);
	if (ndiv2 <= 0)
		ndiv2 = 1;

	get_vars(fp_spc, "total_earthquakes ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &total_earthquakes);

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
	for (int k = 1; k < nxc; k++) {
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
	for (int k = 1; k < nyc; k++) {
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
		for (int k = 1; k < nzc; k++) {
			ib = ie;
			get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
			sscanf(parval, "%d", &igridz[k]);
		}
	}
	fclose(fp_spc);
//-----end of optional parameters

	nx = 1;
	ny = 1;
	nz = 1;

	for (int i = 1; i < nxc; i++) {
		nx = nx + igridx[i - 1];
	}

	for (int i = 1; i < nyc; i++) {
		ny = ny + igridy[i - 1];
	}

	for (int i = 1; i < nzc; i++) {
		nz = nz + igridz[i - 1];
	}
	if(DEBUG_PRINT)
		printf(" Fine grid dimension (nx=%d, ny=%d, nz=%d)\n", nx, ny, nz);

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

// c---If dq and df have not been specified, { make them so that the
// c   interval at the surface is equal to h
// c   First Convert geographic latitude to geocentric colatitude
	if(DEBUG_PRINT)
		printf(" Origin:  x0=%.14lf\t\ty0=%.14lf\t\tz0=%.14lf\n", x0, y[0], z0);
	x00 = x0;
	y00 = y[0];
	y[0] *= degrad;
//   y[0] = hpi - glat(y[0]);
	double z0r;
	y[0] = hpi - glath(y[0], z0, &z0r);
	x0 *= degrad;
	dq *= degrad;
	df *= degrad;
	if (dq == 0)
		dq = h / rearth;
	if (df == 0)
		df = fabs(h / (rearth * sin(y[0])));
	if(DEBUG_PRINT) {
		printf(" Origin:  %.14lf\t%.14lf\t%.14lf\n", x0, y[0] / degrad, z0r);
		printf(" Radial Spacing: %lf\n", h);
		printf(" Latitude Spacing: %.17E\n", dq);
		printf(" Longitude Spacing: %.17E\n", df);
	}
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
	fprintf(fp_log, "          Parameters Set For This Run of Sphfdloc.f\n");
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
		fprintf(fp_log, "\n");

		dtoa(tmp, x00, 18);
		fprintf(fp_log, "  Cartesian X origin (x0):   %s     \n", tmp);
		dtoa(tmp, y00, 18);
		fprintf(fp_log, "  Cartesian Y origin (y0):   %s     \n", tmp);
		dtoa(tmp, z0, 18);
		fprintf(fp_log, "  Cartesian Z origin (z0):   %s     \n", tmp);
		dtoa(tmp, h, 18);
		fprintf(fp_log, "  Fine Radial Spacing     :   %s km  \n", tmp);
		dtoa(tmp, dq / degrad, 18);
		fprintf(fp_log, "  Fine Radial Spacing     :   %s degrees \n", tmp);
		dtoa(tmp, df / degrad, 18);
		fprintf(fp_log, "  Fine Radial Spacing     :   %s degrees \n", tmp);
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

	fprintf(fp_log, "  1st stage grid search division: %12d\n", ndiv);
	fprintf(fp_log, "  2nd stage grid search division: %12d\n", ndiv2);
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Input file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Local Earthquake data file: %s\n", leqsfil);
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Output file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Statistics Summary file: %s\n", fsumfil);
	fprintf(fp_log, " Outlier file: %s\n", outlfil);
	fprintf(fp_log, " New Header file: %s\n", fhedfil);
	fprintf(fp_log, " New Data file: %s\n", fdatfil);
	fprintf(fp_log, " \n");

	fprintf(fp_log, " Travel Time Table Directory: %s\n", timedir);
	fprintf(fp_log, " \n");

	if(DEBUG_PRINT) {
		printf(" Origin:  x0=%.14lf y0=%.14lf z0=%.14lf\n", x0, y[0], z0r);
		printf(" Spacing:  h=%lf\n", h);
		printf(" nxc=%d nyc=%d nzc=%d\n", nxc, nyc, nzc);
	}
	int nxyc = nxc * nyc;
	int nxyzc = nxyc * nzc;

	int nxy = nx * ny;
	int nxyz = nxy * nz;


	float xmax = x0 + (nx - 1) * df;
	float ymax = y[0] + (ny - 1) * dq;
	float zmax = z0 + (nz - 1) * h;
	fprintf(fp_log, " Total Number of fine grid nodes:%13d\n", nxyz);
	fprintf(fp_log, " Total Number of coarse grid nodes:%13d\n", nxyzc);
	fprintf(fp_log, "\n");

	fprintf(fp_log, " X Max:%22.14lf\n", xmax/degrad);
	fprintf(fp_log, " Y Max:%22.14lf\n", ymax/degrad);
	fprintf(fp_log, " Z Max:%22.14lf\n", zmax);
	fprintf(fp_log, "\n");

//---write out some reminders so we know what's going on:
	fprintf(fp_log, " Current settings: \n");
	if (iread == 1) {
		fprintf(fp_log,
				" Test of data read only; no probability execution (iread = 1) \n");
	} else {
		fprintf(fp_log, " Data read as in normal execution (iread = 0)\n");
	}
	if (ivs == 1) {
		fprintf(fp_log, " Ts calculated explicitly (ivs = 1).\n");
	} else {
		fprintf(fp_log, " Ts calculated as Tp*vpvs (ivs = 0). vpvs = %lf\n",
				vpvs);
	}
	fprintf(fp_log, " Default value for Vp/Vs = %lf\n", vpvs);
	fprintf(fp_log, " Number of Phases Threshold    : %d\n", nthres);
	fprintf(fp_log, " Absolute Residual Threshold   : %lf\n", resthres);
	fprintf(fp_log, " Percentage Residual Threshold : %lf\n", resthrep);
	fprintf(fp_log, " Standard Deviation Threshold  : %lf\n", stdmax);
	fprintf(fp_log, " Beginning Depth grid          : %d\n", kmin + 1);
	fprintf(fp_log, "\n");
	fprintf(fp_log,
			"*************************************************************** \n");
	fprintf(fp_log, "\n");

//---START LOOP OVER EVENTS
//   read in event header
	if(total_earthquakes == 0) {
		total_earthquakes = earthquake_file_delimiter(leqsfil, eqkdir);
	}
	char timefiles[maxsta][MAXSTRLEN + 1];
	int timefile_counts = read_timefiles(iread, nxyz, timefiles);
	if (timefile_counts < 0) {
		printf("file can not open\n");
		assert(0);
	}
	//qsort(timefiles, timefile_counts, sizeof(timefiles[0]), (int (*)(const void*, const void*))strcmp);
	/*
	char *str_fhd[total_earthquakes], *str_data[total_earthquakes];
	char *str_sum[total_earthquakes], *str_out[total_earthquakes];
	for(int i=0;i<total_earthquakes;i++) {
		str_fhd[i]=(char*)malloc(sizeof(char)*100);
		str_out[i]=(char*)malloc(sizeof(char)*1000);
		str_data[i]=(char*)malloc(sizeof(char)*10000);
		str_sum[i]=(char*)malloc(sizeof(char)*20000);
		memset(str_fhd[i], 0, sizeof(str_fhd[i]));
		memset(str_data[i], 0, sizeof(str_data[i]));
		memset(str_sum[i], 0, sizeof(str_sum[i]));
		memset(str_out[i], 0, sizeof(str_out[i]));
	}
	*/

//#pragma omp parallel for
	for (int nev = 0; nev < total_earthquakes; nev++) {
		char filename[100];
		sprintf(filename, "%s/%d.eqk", eqkdir, nev + 1);
		FILE *fp_leq = fopen(filename, "r");
		if (!fp_leq) {
			printf("file not found: %s\n", filename);
			assert(0);
		}
		char str_inp[MAXSTRLEN];
		if (fgets(str_inp, sizeof(str_inp), fp_leq) == 0) {
			if (feof(fp_leq)) {
				printf("eof on fp_leq(%s)\n", filename);
				assert(0);
			} else {
				printf("error on reading fp_leq(%s)\n", filename);
				assert(0);
			}
		}
		if(strlen(str_inp) > sizeof(str_inp)) {
			printf("str_inp is too long: %s\n", str_inp);
			assert(0);
		}
		int iyr, jday, ihr, imn;
		float sec;
		float xlat, xlon, dep;
		char evid[10 + 1];
// c---START LOOP OVER EVENTS
// c	read in event header
		sscanf(str_inp, "%d %d %d %d %f %f %f %f %s\n", &iyr, &jday, &ihr,
				&imn, &sec, &xlat, &xlon, &dep, evid);
//	c---convert to epochal time
		double dsec = sec;
		double dpot;	
		htoe(iyr, jday, ihr, imn, dsec, &dpot);

		int isgood[maxobs], indsta[maxobs];
		float obstime[maxobs], res[maxobs], resmin[maxobs], wtmin[maxobs];
		float wtsave[maxobs], pwt[maxobs], rwts[maxobs], tps[maxobs];

		memset(isgood, 0, sizeof(isgood));
		memset(indsta, 0, sizeof(indsta));

		char sta[maxobs][MAXSTRLEN + 1];
		char phs[maxobs];
		char usemark;
		int nsta = 0;
// c	read in the times and stations
		for(nsta = 0; fgets(str_inp, sizeof(str_inp), fp_leq); nsta++) {
			if (nsta >= maxobs) {
				printf("Error:  too many observations! nsta=%d maxobs=%d\n",
						nsta, maxobs);
				assert(0);
			}
			
			trim(str_inp);
			int len = strlen(str_inp);
			if (str_inp[0] == '\n') {
				break;
			}
			if (len > MAXSTRLEN) {
				printf("input length is too large. len=%d str_inp=%s\n", len,
						str_inp);
				assert(0);
			}

			char str_tmp[100];
			sscanf(str_inp, "%s %d %d %d %d %f %99[^\n]\n", sta[nsta], &iyr,
					&jday, &ihr, &imn, &sec, str_tmp);

			isgood[nsta] = 1;
			if (str_tmp[0] == '*') {
				isgood[nsta] = 0;
				usemark = '*';
				str_tmp[0] = ' ';
				trim(str_tmp);
			}
			sscanf(str_tmp, "%c %f", &phs[nsta], &rwts[nsta]);
// c---convert from "0 1" to "P S" if necessary
			if (phs[nsta] == '0') {
				phs[nsta] = 'P';
			} else if (phs[nsta] == '1') {
				phs[nsta] = 'S';
			}
			pwt[nsta] = 1.f / (rwts[nsta] * rwts[nsta]);
			dsec = sec;

			double tarr;
			htoe(iyr, jday, ihr, imn, dsec, &tarr);
			obstime[nsta] = tarr - dpot;
			char key[MAXSTRLEN];
			strcpy(key, sta[nsta]);
			if (phs[nsta] == 'P') {
				strcat(key, ".ptimes");
			} else if (phs[nsta] == 'S') {
				strcat(key, ".stimes");
			} else {
				assert(0);
			}
			//char *word = (char*)bsearch(key, timefiles, timefile_counts, MAXSTRLEN, (int(*)(const void*, const void*)) strcmp);

			indsta[nsta] = -1;
			for (int ind = 0; ind < timefile_counts; ind++) {
				if (strcmp(timefiles[ind], key) == 0) {
					indsta[nsta] = ind;
					break;
				}
			}
			if (indsta[nsta] == -1) {
				printf("file: %s not found!\n", key);
				assert(0);
			}
		}
		fclose(fp_leq);
		if(DEBUG_PRINT)
			printf(" %d times read in for event %d\n", nsta, nev + 1);

		int ngood = nsta;
		if(DEBUG_PRINT)
			printf(" Working on event %d\n", nev + 1);

//---see if the travel time tables for these stations have been read in.  If
//   not, read them in.
		if (nsta > maxsta) {
			printf("Error: too many station.\n");
			assert(0);
		}

		char pha[maxsta];
		char stn[maxsta][MAXSTRLEN];
		int ntread = 0;
		for (int i = 0; i < nsta; i++) {
			if(DEBUG_PRINT)
				printf("%4d%4d %4s %c\n", i + 1, ntread, sta[i], phs[i]);
			int j;
			if (ntread > 0) {
				for (j = 0; j < ntread; j++) {
					if (strcmp(sta[i], stn[j]) == 0) {
						if (ivs == 0 || phs[i] == pha[j]) {
							continue;
						}
					}
				}
			}
			int iuse = ntread;
			ntread++;
			if (ntread > maxsta) {
				if(DEBUG_PRINT)
					printf(" Limit reached...looking for one we do not need now.\n");
				assert(0);
				for (j = 0; j < ntread; j++) {
					int flag = 0;
					for (int k = 0; k < nsta; k++) {
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
				if(DEBUG_PRINT)
					printf("Sorry, no more space available...\n");
				assert(0);
				a21: ntread--;
				iuse = j;
				if(DEBUG_PRINT)
					printf("Overwriting %s  %c with  %s %c\n", stn[j], pha[j],
						sta[i], phs[i]);
			}
			strcpy(stn[iuse], sta[i]);
			pha[iuse] = phs[i];
		}
		if(DEBUG_PRINT)
			printf(" Current ntread = %d\n", ntread);

		//----test to see if all data can be read in correctly
		if (iread == 1)
			continue;

		//   Start the PDF calculation.  Loop over all grid points.
		int iretry = 0;
		int maxretry = 1;

		a33: 
		if(DEBUG_PRINT)
			printf(" Starting Coarse Grid Loop ... \n");
		//---NB:  we start k at an nz corresponding to a minimum depth.
		//   h = 10, z0 = -25, so starting at k = 4 makes min depth at this stage 5.0.
		//   this could be anywhere from -5 to 5 in next stage.
		int nxx = 0, nyy = 0, nzz = 0;
		float stdmin = DBL_MAX;
		float avrmin = 0.;
		for (int k = kmin; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				for (int i = 0; i < nx; i++) {
					int n = nxy * k + nx * j + i;
					float avres = 0.0, rsq = 0.0, facs = 0.0;
					for (int is = 0; is < nsta; is++) {
						double tp = t[indsta[is]][n];
						if (ivs == 0 && phs[is] == 'S') {
							tp *= vpvs;
						}
						res[is] = obstime[is] - tp;
						float wt = pwt[is];
						wtsave[is] = wt;
						if (isgood[is]) {
							float res1 = res[is] * wt;
							float res2 = res1 * res[is];
							avres += res1;
							rsq += res2;
							facs += wt;
						}
					}
					//    c---calculate variance of data
					//    c         avwt = facs/nsta
					//    c         deb = (rsq - avres*avres/facs)/((nsta - 4)*avwt)
					float avwt = facs / ngood;
					float deb = (rsq - avres * avres / facs)
							/ ((ngood - 4) * avwt);
					float dsum = fabs(deb);
					//    c------ find total sigma and the most likely hypocenter with min std
					float std = sqrtf(dsum);
					if (std < stdmin) {
						nxx = i;
						nyy = j;
						nzz = k;
						stdmin = std;
						avrmin = avres / facs;
						for (int i2 = 0; i2 < nsta; i2++) {
							resmin[i2] = res[i2];
							wtmin[i2] = wtsave[i2];
							tps[i2] = obstime[i2] - res[i2];
						}
					}
				}
			}
		}

		double exm = x0 + nxx * df;
		double eym = y[0] + nyy * dq;
		double ezm = z0 + nzz * h;

		xlon = exm / degrad;
		//   xlat = glatinv(hpi-eym)/degrad
		double ezmr = 0;
		xlat = glathinv(hpi - eym, z0r - (ezm - z0), &ezmr) / degrad;
		if(DEBUG_PRINT) {
			printf("\n");
			printf("\n");
			printf(" Coarse grid finished ... \n");
			printf(" Info Density Maximum for event %12d\n", nev + 1);
			printf(" Minimum Standard Deviation   : %12.8lf\n", stdmin);
			printf(" Optimal Grid point (nx,ny,nz): %12d %11d %11d\n", nxx + 1,
					nyy + 1, nzz + 1);
			printf(" Average residual             : %10.8lf\n", avrmin);
			printf(" Hypocenter (lat, lon ,z )    :   %12.3f  %12.3f  %12.3f\n",
					xlat, xlon, ezmr);
			printf("\n");
			printf(" Station statistics at this point: \n");
			printf(
					"Code  Ph  Demeaned    True   Expected  Used   Tobs     Tcalc     To-Tc \n");
			printf(
					"          Residual   Weight   Uncert.  Flag    Abs      Abs     \n");
			printf("\n");
			for (int i = 0; i < nsta; i++) {
				printf("%4s   %c %9.4lf %9.3E %9.4f %3d %9.4f %9.4f %9.4f\n",
						sta[i], phs[i], resmin[i] - avrmin, wtmin[i],
						1.f / sqrt(wtmin[i]), isgood[i], obstime[i], tps[i],
						obstime[i] - tps[i]);
			}
			printf("\n");
			printf("\n");
		}
		//----Now refine location to subgrid spacing
		//   ndiv = 20
		int nstep1 = ndiv + 1;
		int nstep2 = nstep1 + ndiv;
		double dh = h / ndiv;
		double ddq = dq / ndiv;
		double ddf = df / ndiv;

		int iz0 = nzz - 1;
		int nstepz = nstep2;
		//   if (iz0.le.0) {
		if (iz0 < kmin) {
			iz0 = nzz;
			nstepz = nstep1;
		}
		if (nzz == nz - 1) {
			nstepz = nstep1 - 1;
		}

		int iy0 = nyy - 1;
		int nstepy = nstep2;
		if (iy0 < 0) {
			iy0 = nyy;
			nstepy = nstep1;
		}
		if (nyy == ny - 1) {
			nstepy = nstep1 - 1;
		}

		int ix0 = nxx - 1;
		int nstepx = nstep2;
		if (ix0 < 0) {
			ix0 = nxx;
			nstepx = nstep1;
		}
		if (nxx == nx - 1) {
			nstepx = nstep1 - 1;
		}

		double z0i = z0 + h * iz0;
		double y0i = y[0] + dq * iy0;
		double x0i = x0 + df * ix0;
		if(DEBUG_PRINT)
			printf("ix0=%d iy0=%d iz0=%d\n", ix0, iy0, iz0);
		int count = 0;
		for (int k = 0; k < nstepz; k++) {
			double zp = z0i + dh * k;
			for (int j = 0; j < nstepy; j++) {
				double yp = y0i + ddq * j;
				for (int i = 0; i < nstepx; i++) {
					float avres = 0.0, rsq = 0.0, facs = 0.0;
					double xp = x0i + ddf * i;
					count ++;
					for (int is = 0; is < nsta; is++) {
						double tp = 0;
						find_time(xp, yp, zp, &tp, is, indsta);
						//printf("%d %lf\n", count, tp);
						if (ivs == 0 && phs[is] == 'S') {
							tp *= vpvs;
						}
						res[is] = obstime[is] - tp;
						float wt = pwt[is];
						wtsave[is] = wt;
						if (isgood[is]) {
							float res1 = res[is] * wt;
							float res2 = res1 * res[is];
							avres += res1;
							rsq += res2;
							facs += wt;
						}
					}
					//    c---calculate variance of data
					//    c         avwt = facs/nsta
					//    c         deb = (rsq - avres*avres/facs)/((nsta - 4)*avwt)
					float avwt = facs / ngood;
					float deb = (rsq - avres * avres / facs)
							/ ((ngood - 4) * avwt);
					float dsum = fabs(deb);
					//------ find total sigma and the most likely hypocenter with min std
					float std = sqrtf(dsum);
					if (std < stdmin) {
						nxx = i;
						nyy = j;
						nzz = k;
						exm = xp;
						eym = yp;
						ezm = zp;
						stdmin = std;
						avrmin = avres / facs;
						for (int i2 = 0; i2 < nsta; i2++) {
							resmin[i2] = res[i2];
							tps[i2] = obstime[i2] - res[i2];
							wtmin[i2] = wtsave[i2];
						}
					}
				}
			}
		}
		if(DEBUG_PRINT) {
			printf("\n");
			printf("\n");
			printf(" First Order Refinement finished ... \n");
			printf(" Info Density Maximum for event   %12d\n", nev + 1);
			printf(" Minimum Standard Deviation   : %12.9lf\n", stdmin);
			printf(" Optimal Grid point (nx,ny,nz): %12d %11d %11d\n", nxx + 1,
					nyy + 1, nzz + 1);
			printf(" Average residual             : %10.8lf\n", avrmin);
			printf(" Hypocenter (lat, lon ,z )    :   %12.3f  %12.3f  %12.3f\n",
					exm, eym, ezm);
			printf("\n");
			printf(" Station statistics at this point: \n");
			printf(
					"Code  Ph  Demeaned    True   Expected  Used   Tobs     Tcalc     To-Tc \n");
			printf(
					"          Residual   Weight   Uncert.  Flag    Abs      Abs     \n");
			for (int i = 0; i < nsta; i++) {
				printf("%4s   %c %9.4lf %9.3E %9.4f %3d %9.4f %9.4f %9.4f\n",
						sta[i], phs[i], resmin[i] - avrmin, wtmin[i],
						1. / sqrtf(wtmin[i]), isgood[i], obstime[i], tps[i],
						obstime[i] - tps[i]);
			}
			printf("\n");
			printf("\n");
		}
		//----Second order refinement
		//   ndiv2 = 20
		nstep1 = ndiv2 + 1;
		nstep2 = nstep1 + ndiv2;
		double ddh = dh / ndiv2;
		double dddq = ddq / ndiv2;
		double dddf = ddf / ndiv2;

		z0i = ezm - dh;
		nstepz = nstep2;
		//   if (z0i.le.z0) {
		if (z0i <= z0 + kmin * h) {
			z0i = ezm;
			nstepz = nstep1;
		}
		if (ezm + dh > zmax)
			nstepz = nstep1 - 1;

		y0i = eym - ddq;
		nstepy = nstep2;
		if (y0i <= y[0]) {
			y0i = eym;
			nstepy = nstep1;
		}
		if (eym + ddq > ymax)
			nstepy = nstep1 - 1;

		x0i = exm - ddf;
		nstepx = nstep2;
		if (x0i <= x0) {
			x0i = exm;
			nstepx = nstep1;
		}

		if (exm + ddf > xmax)
			nstepx = nstep1 - 1;
		for (int k = 0; k < nstepz; k++) {
			double zp = z0i + ddh * k;
			for (int j = 0; j < nstepy; j++) {
				double yp = y0i + dddq * j;
				for (int i = 0; i < nstepx; i++) {
					float avres = 0.0, rsq = 0.0, facs = 0.0;
					double xp = x0i + dddf * i;
					for (int is = 0; is < nsta; is++) {
						double tp = 0;
						find_time(xp, yp, zp, &tp, is, indsta);
						if (ivs == 0 && phs[is] == 'S') {
							tp *= vpvs;
						}
						res[is] = obstime[is] - tp;
						float wt = pwt[is];
						wtsave[is] = wt;
						if (isgood[is]) {
							float res1 = res[is] * wt;
							float res2 = res1 * res[is];
							avres += res1;
							rsq += res2;
							facs += wt;
						}
					}
					//    c---calculate variance of data
					//    c         avwt = facs/nsta
					//    c         deb = (rsq - avres*avres/facs)/((nsta - 4)*avwt)
					float avwt = facs / ngood;
					float deb = (rsq - avres * avres / facs)
							/ ((ngood - 4) * avwt);
					float dsum = fabs(deb);
					//------ find total sigma and the most likely hypocenter with min std
					float std = sqrtf(dsum);
					if (std < stdmin) {
						nxx = i;
						nyy = j;
						nzz = k;
						exm = xp;
						eym = yp;
						ezm = zp;
						stdmin = std;
						avrmin = avres / facs;
						for (int i2 = 0; i2 < nsta; i2++) {
							resmin[i2] = res[i2];
							tps[i2] = obstime[i2] - res[i2];
							wtmin[i2] = wtsave[i2];
						}
					}
				}
			}
		}

		xlon = exm / degrad;
		//   xlat = glatinv(hpi-eym)/degrad
		xlat = glathinv(hpi - eym, z0r - (ezm - z0), &ezmr) / degrad;
		if(DEBUG_PRINT) {
			printf("xlon=%lf xlat=%lf\n", xlon, xlat);
			printf("\n");
			printf("\n");
			printf(" Second Order Refinement finished ... \n");
			printf(" Info Density Maximum for event   %12d\n", nev + 1);
			printf(" Minimum Standard Deviation   : %12.9lf\n", stdmin);
			printf(" Optimal Grid point (nx,ny,nz): %12d %11d %11d\n", nxx + 1,
					nyy + 1, nzz + 1);
			printf(" Average residual             : %10.8lf\n", avrmin);
			printf(" Hypocenter (lat, lon ,z )    :   %12.3f  %12.3f  %12.3f\n",
					xlat, xlon, ezmr);
			printf("\n");
			printf(" Station statistics at this point: \n");
			printf(
					"Code  Ph  Demeaned    True   Expected  Used   Tobs     Tcalc     To-Tc \n");
			printf(
					"          Residual   Weight   Uncert.  Flag    Abs      Abs     \n");
			for (int i = 0; i < nsta; i++) {
				printf("%4s   %c %9.4lf %9.3E %9.4f %3d %9.4f %9.4f %9.4f\n",
						sta[i], phs[i], resmin[i] - avrmin, wtmin[i],
						1. / sqrtf(wtmin[i]), isgood[i], obstime[i], tps[i],
						obstime[i] - tps[i]);
			}
			printf("\n");
			printf("\n");
		}
		//----check for outliers and acceptable locations
		ngood = 0;
		int iredo = 0;
		for (int i = 0; i < nsta; i++) {
			if (isgood[i]) {
				float absres = fabs(resmin[i] - avrmin);
				float percres = 100.0f * absres / tps[i];
				if (absres <= resthres || percres <= resthrep)
					ngood++;
				else {
					isgood[i] = 0;
					iredo = 1;
				}
			}
		}
		int len_str_out = 0;
		if (ngood >= nthres) {
			//---redo the location if outliers have been removed from a well recorded event
			if (iredo == 1) {
				if(DEBUG_PRINT)
					printf(" Outliers found .. redoing location \n");
				goto a33;
			}
			//---in some cases, a revised location can bring back data that originally was considered outlier
			//       so we allow the program to go back and use this recovered data
			if (iretry < maxretry) {
				iretry++;
				iredo = 0;
				for (int i = 0; i < nsta; i++) {
					if (isgood[i] == 0) {
						float absres = fabs(resmin[i] - avrmin);
						float percres = 100.0f * absres / tps[i];
						if (absres <= resthres || percres <= resthrep) {
							ngood++;
							isgood[i] = 1;
							iredo = 1;
						}
					}
				}
				if (iredo == 1) {
					if(DEBUG_PRINT)
						printf(" Recovered data .. redoing location \n");
					goto a33;
				}
			}
			if (stdmin <= stdmax) {
				xlon = exm / degrad;
				//       xlat = glatinv[hpi-eym]/degrad;
				xlat = glathinv(hpi - eym, z0r - (ezm - z0), &ezmr) / degrad;
				double ot = dpot + avrmin;
				etoh(ot, &iyr, &jday, &ihr, &imn, &dsec);
				sec = dsec;
				char tmp[100];
				sprintf(tmp,
						"%4d %3d %2d %2d %8.4f %9.5f %10.5f %8.4lf %12s %13.3lf\n",
						iyr, jday, ihr, imn, sec, xlat, xlon, ezmr, evid,
						stdmin);
				int len_str_data = 0;
				strapp(str_data[nev], &len_str_data, tmp);
				strcpy(str_fhd[nev], tmp);
				
				sprintf(tmp,
						"%4d %3d %2d %2d %8.4f %9.5f %10.5f %8.4lf %12s %13.3f\n",
						iyr, jday, ihr, imn, sec, exm, eym, ezmr, evid, stdmin);

				if (ngood < nsta) {
					strapp(str_out[nev], &len_str_out, tmp);
				}
				if(DEBUG_PRINT) {
					printf("\n");
					printf(tmp);

					printf(" Event ID: %s Finished\n", evid);
					printf("\n");
				}
				int len_str_sum = 0;
				strapp(str_sum[nev], &len_str_sum, "\n\n");
				sprintf(tmp, "  Minimum Standard Deviation   : %12.8lf\n",
						stdmin);
				strapp(str_sum[nev], &len_str_sum, tmp);
				sprintf(tmp, "  Info Density Maximum for event %12d\n", nev);
				strapp(str_sum[nev], &len_str_sum, tmp);
				sprintf(tmp,
						"%4d %3d %2d %2d %8.4f  %8.5f %10.5f %8.4lf %12s %13.3lf\n",
						iyr, jday, ihr, imn, sec, exm, eym, ezmr, evid, stdmin);
				strapp(str_sum[nev], &len_str_sum, tmp);
				sprintf(tmp, " Event ID: %7s\n", evid);
				strapp(str_sum[nev], &len_str_sum, tmp);
				sprintf(tmp, "  Minimum Standard Deviation   : %12.8lf\n",
						stdmin);
				strapp(str_sum[nev], &len_str_sum, tmp);
				strcpy(tmp,
						"  Maximum Sigma                :   2.51547135E+15\n");
				strapp(str_sum[nev], &len_str_sum, tmp);
				strcpy(tmp,
						"  Total Sigma                  :    0.0000000    \n");
				strapp(str_sum[nev], &len_str_sum, tmp);
				sprintf(tmp,
						"  Optimal Grid point (nx,ny,nz): %12d %11d %11d\n",
						nxx + 1, nyy + 1, nzz + 1);
				strapp(str_sum[nev], &len_str_sum, tmp);
				sprintf(tmp, "  Average residual             :   %10.8lf\n",
						avrmin);
				strapp(str_sum[nev], &len_str_sum, tmp);
				sprintf(tmp,
						"  Hypocenter (lat, lon ,z )    :   %12.3f  %12.3f  %12.3f\n\n",
						xlat, xlon, ezmr);
				strapp(str_sum[nev], &len_str_sum, tmp);
				strcpy(tmp, "  Station statistics at this point: \n");
				strapp(str_sum[nev], &len_str_sum, tmp);
				strcpy(tmp,
						" Code  Ph  Demeaned    True   Expected  Used   Tobs     Tcalc     To-Tc \n");
				strapp(str_sum[nev], &len_str_sum, tmp);
				strcpy(tmp,
						"           Residual   Weight   Uncert.  Flag    Abs      Abs     \n\n");
				strapp(str_sum[nev], &len_str_sum, tmp);
				for (int i = 0; i < nsta; i++) {
					sprintf(tmp,
							"%4s   %c %10.4lf %9.3E %9.4f %3d %9.4f %9.4f %9.4f\n",
							sta[i], phs[i], resmin[i] - avrmin, wtmin[i],
							1. / sqrtf(wtmin[i]), isgood[i], obstime[i], tps[i],
							obstime[i] - tps[i]);
					strapp(str_sum[nev], &len_str_sum, tmp);
				}
				strapp(str_sum[nev], &len_str_sum, "\n");
				for (int j = 0; j < nsta; j++) {
					double tarr = obstime[j] + dpot;
					etoh(tarr, &iyr, &jday, &ihr, &imn, &dsec);
					if (isgood[j]) {
						sprintf(tmp,
								"%-6s%4d %3d %2d %2d %7.3lf %c %15.3f %10.3lf\n",
								sta[j], iyr, jday, ihr, imn, dsec, phs[j],
								rwts[j], resmin[j] - avrmin);
					} else {
						sprintf(tmp,
								"%-6s%4d %3d %2d %2d %7.3lf*%c %15.3f %10.3lf\n",
								sta[j], iyr, jday, ihr, imn, dsec, phs[j],
								rwts[j], resmin[j] - avrmin);
						strapp(str_out[nev], &len_str_out, tmp);
					}
					strapp(str_data[nev], &len_str_data, tmp);
				}
				strapp(str_data[nev], &len_str_data, " \n");
			}
		} else {
			double ot = dpot + avrmin;
			etoh(ot, &iyr, &jday, &ihr, &imn, &dsec);
			sec = dsec;
			char tmp[200];
			if (DBL_MAX - stdmin < 1) {
				sprintf(tmp,
						"%4d %3d %2d %2d %8.4f %9.5f %10.5lf %8.4lf %12s    **********\n",
						iyr, jday, ihr, imn, sec, xlat, xlon, ezmr, evid);
			} else {
				sprintf(tmp,
						"%4d %3d %2d %2d %8.4f %9.5f %10.5lf %8.4lf %12s %13.3lf\n",
						iyr, jday, ihr, imn, sec, xlat, xlon, ezmr, evid,
						stdmin);
			}
			strapp(str_out[nev], &len_str_out, tmp);
			for (int j = 0; j < nsta; j++) {
				float tarr = obstime[j] + dpot;
				etoh(tarr, &iyr, &jday, &ihr, &imn, &dsec);
				sprintf(tmp, "%s %4d %3d %2d %2d %7.3lf*%c %15.3f %10.3lf\n",
						sta[j], iyr, jday, ihr, imn, dsec, phs[j], rwts[j],
						resmin[j] - avrmin);
				strapp(str_out[nev], &len_str_out, tmp);
			}
			strapp(str_out[nev], &len_str_out, "\n");
		}
	}
//---END OF LOOP OVER EQS
	FILE *fp_fhd = fopen(fhedfil, "w");
	FILE *fp_fdt = fopen(fdatfil, "w");
	FILE *fp_sum = fopen(fsumfil, "w");
	FILE *fp_out = fopen(outlfil, "w");

	if (!fp_fhd) {
		printf("fp_fhd: open '%s' error\n", fhedfil);
		assert(0);
	}
	if (!fp_fdt) {
		printf("fp_fdt: open '%s' error\n", fdatfil);
		assert(0);
	}
	if (!fp_sum) {
		printf("fp_sum: open '%s' error\n", fsumfil);
		assert(0);
	}
	if (!fp_out) {
		printf("fp_out: open '%s' error\n", outlfil);
		assert(0);
	}

	for (int i = 0; i < total_earthquakes; i++) {
		if (fprintf(fp_fhd, "%s", str_fhd[i]) < 0) {
			printf("write fp_fhd(%s) error\n", fhedfil);
			assert(0);
		}
		if (fprintf(fp_fdt, "%s", str_data[i]) < 0) {
			printf("write fp_fdt(%s) error\n", fdatfil);
			assert(0);
		}
		if (fprintf(fp_sum, "%s", str_sum[i]) < 0) {
			printf("write fp_fdt(%s) error\n", fsumfil);
			assert(0);
		}
		if (fprintf(fp_out, "%s", str_out[i]) < 0) {
			printf("write fp_fdt(%s) error\n", outlfil);
			assert(0);
		}
	}

	fclose(fp_fhd);
	fclose(fp_fdt);
	fclose(fp_sum);
	fclose(fp_out);
	return 0;
}

void find_time(double x, double yy, double z, double *tp, int is, int *indsta) {
	int i = ((int) ((x - x0) / df));
	int j = ((int) ((yy - y[0]) / dq));
	int k = ((int) ((z - z0) / h));
	//printf("i,j,k = %d %d %d \n",i,j,k);
	if (i < 0)
		i = 0;
	if (j < 0)
		j = 0;
	if (k < 0)
		k = 0;

	if (i > nx - 2)
		i = nx - 2;
	if (j > ny - 2)
		j = ny - 2;
	if (k > nz - 2)
		k = nz - 2;

	double xi = x0 + df * i;
	double yj = y[0] + dq * j;
	double zk = z0 + h * k;

	int nk = nx * ny * k;
	int nj = nx * j;
	int nk2 = nx * ny * ( k + 1);
	int nj2 = nx * (j + 1);
	double fx = (x - xi) / df;
	double fy = (yy - yj) / dq;
	double fz = (z - zk) / h;

	int ist = indsta[is];
	*tp = (1.f - fx) * (1.f - fy) * (1.f - fz) * t[ist][nk + nj + i];
	*tp = *tp + fx * (1.f - fy) * (1.f - fz) * t[ist][nk + nj + i + 1];
	*tp = *tp + (1.f - fx) * fy * (1.f - fz) * t[ist][nk + nj2 + i];
	*tp = *tp + (1.f - fx) * (1.f - fy) * fz * t[ist][nk2 + nj + i];
	*tp = *tp + fx * fy * (1.f - fz) * t[ist][nk + nj2 + i + 1];
	*tp = *tp + fx * (1.f - fy) * fz * t[ist][nk2 + nj + i + 1];
	*tp = *tp + (1.f - fx) * fy * fz * t[ist][nk2 + nj2 + i];
	*tp = *tp + fx * fy * fz * t[ist][nk2 + nj2 + i + 1];
}

void read_station_set(int *nsta, int *iyr, int *jday, int *ihr, int *imn,
		float *sec, int *isgood, char *phs, float *rwts,
		char sta[maxobs][MAXSTRLEN + 1], char *usemark, FILE *fp_leq) {
	while (1) {
		char str_inp[100];
		if (fgets(str_inp, sizeof(str_inp), fp_leq) == 0) {
			printf("reading fp_leq error.\n");
			assert(0);
		}
		trim(str_inp);
		int len = strlen(str_inp);
		if (str_inp[0] == '\n') {
			break;
		}
		printf("str_inp=%s\n", str_inp);
		if (str_inp[len - 1] != '\n') {
			printf("input length is too large. len=%d str_inp=%s\n", len,
					str_inp);
			assert(0);
		}

		char str_tmp[100];
		sscanf(str_inp, "%s %d %d %d %d %f %99[^\n]\n", sta[*nsta],
				&iyr[*nsta], &jday[*nsta], &ihr[*nsta], &imn[*nsta],
				&sec[*nsta], str_tmp);

		isgood[*nsta] = 1;
		if (str_tmp[0] == '*') {
			isgood[*nsta] = 0;
			*usemark = '*';
			str_tmp[0] = ' ';
			trim(str_tmp);
		}
		sscanf(str_tmp, "%c %f", &phs[*nsta], &rwts[*nsta]);
		*nsta = *nsta + 1;
		if (*nsta > maxobs) {
			printf("Error:  too many observations! nsta=%d maxobs=%d\n", *nsta,
					maxobs);
			assert(0);
		}
	}
}

int read_timefiles(int iread, int nxyz, char timefiles[maxsta][MAXSTRLEN + 1]) {
	DIR *d = opendir(timedir);
	if (!d) {
		printf("opendir error: %s\n", timedir);
		return -1;
	}

	struct dirent *dir;
	int i = -1;
	for (i = 0; ((dir = readdir(d)) != NULL); i++) {
		char *filename = dir->d_name;
		if (filename == NULL) {
			printf("FileName is NULL\n");
			assert(0);
		}
		if ((strstr(filename, ".ptimes") == NULL)
				&& (strstr(filename, ".stimes") == NULL)) {
			i--;
			continue;
		}

		//---see if the travel time tables for these stations have been read in.  If not, read them in.
		if (i >= maxsta) {
			printf("Error: too many station.\n");
			assert(0);
		}

		strcpy(timefiles[i], filename);
		char filen[80 + 1];
		sprintf(filen, "%s/%s", timedir, filename);
		printf("Reading in %s\n", filen);

// c---see if this the time file has a header on it.
		FILE *fp_tab = fopen(filen, "rb");
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
		sscanf(offset, "%124s", hcomm);

//---verify that this is a valid header
		if (strcmp(head, "HEAD") == 0) {
			if(DEBUG_PRINT) {
				printf(" File has a header...\n");
			}
			if (strcmp(type, "FINE") != 0) {
				if(DEBUG_PRINT)
					printf("WARNING: input mesh does not appear to be FINE: %s\n", type);
			}
			if (iread == 0) {
				fread(t[i], sizeof(t[i][0]), nxyz*4, fp_tab);
			}
		} else {
			rewind(fp_tab);
			if(DEBUG_PRINT) {
				printf(" File has no header...\n");
			}
			if (iread == 0) {
				fread(t[i], sizeof(t[i][0]), nxyz*4, fp_tab);
			}
		}
		fclose(fp_tab);
		if(DEBUG_PRINT)
			printf("...Done.\n");
	}
	return i;
}
