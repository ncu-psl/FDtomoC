#include "common/read_spec.h"
#include "common/gridspec.h"
#include "common/shared_variables.h"

#define MUSTV  4
#define MUSTF  21
#define MAXSTRLEN 132

int nxc, nyc, nzc, nx, ny, nz;
double h, x0, y[1], z0, dq, df, x00, y00;
int igridx[nxcm1], igridy[nycm1], igridz[nzcm1];

double clat, clon, cz;
float az, azmod;
int iflat = 0, isph = 0, vs1d = 1;

//sphfdloc
int iread = 0;
int ivs = 1;
double vpvs = 1.78;
int nthres = 8;
double resthres = .5;
double resthrep = 5.0;
double stdmax = 15.0;
int kmin = 2;
int ndiv = 20;
int ndiv2 = 20;
int ittnum = 1;
int total_earthquakes = 0;
char timedir[60 + 1] = "./ttimes\0";
char eqkdir[60 + 1] = "./eqkdir/\0";

//sphrayderv
float vpvsscale = 0;
int idmean = 0, iray = 0, iraystat = 0, idatout = 1, nomat = 0;
int ivpvs = 0, istacor = 0, idoshot = 0, idotel = 0, kmax;
float resflag = 1.0;

//runlsqr
float damper = 0.001;
int intlims = 0;

//makenewmod
int limitu = 0, mavx = 3, mavy = 3, mavz = 3;
int nsmooth = 2, ipscflg = 0, ido1d = 0;
float dvperc = 0.05, pertscl = 1.0;

//files
char oldvfil[MAXSTRLEN + 1], onedfil[MAXSTRLEN + 1]; // make1d
char tgrdfil[MAXSTRLEN + 1], finevel[MAXSTRLEN + 1]; // c2f
char leqsfil[MAXSTRLEN + 1], fsumfil[MAXSTRLEN + 1], outlfil[MAXSTRLEN + 1],
             fhedfil[MAXSTRLEN + 1], fdatfil[MAXSTRLEN + 1]; //sphfdloc
char stafile[MAXSTRLEN + 1], locdfil[MAXSTRLEN + 1], telrerr[MAXSTRLEN + 1], 
            dtdsfil[MAXSTRLEN + 1], resfile[MAXSTRLEN + 1], hitfile[MAXSTRLEN + 1], 
			dtdhfil[MAXSTRLEN + 1], bookfil[MAXSTRLEN + 1], sclefil[MAXSTRLEN + 1];  //sphrayderv
char nmodfil[MAXSTRLEN + 1], fresfil[MAXSTRLEN + 1]; //runlsqr			
char fmodfil[MAXSTRLEN + 1]; //makenewmod
void read_variables(char *spec_file){
    char *mvals[MUSTV] = { "nxc\0", "nyc\0", "nzc\0", "h\0" };
	char pval[MAXSTRLEN + 1];
    int len, ierr;
    FILE *fp_spc;

    fp_spc = fopen(spec_file, "r");
	if (!fp_spc) {
		printf("(Error in read_spec.c)read fp_spc file error.\n");
		assert(0);
	}
	

//---recover the variables needed to run this program
//
//       nxc, nyc, nzc      coarse dimensions of the fine mesh used in the trt tables
//       h                  fine grid spacing
    int i;
	for (i = 0; i < MUSTV; i++) {
		get_vars(fp_spc, mvals[i], pval, &len, &ierr);
		if (ierr == 1) {
			read_error(mvals[i], "variable", fp_spc);
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
	get_vars(fp_spc, "azmod ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &azmod);
	get_vars(fp_spc, "ittnum ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &ittnum);
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

	//sphfdloc
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
		sscanf(pval, "%lf", &vpvs);
	get_vars(fp_spc, "ivpvs ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%d", &ivpvs);
	}
	if (ivpvs != 0 && ivpvs != 1) {
		ivpvs = 0;
	}

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

	get_vars(fp_spc, "nthres ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &nthres);
	get_vars(fp_spc, "resthres ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &resthres);
	get_vars(fp_spc, "resthrep ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &resthrep);
	get_vars(fp_spc, "stdmax ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &stdmax);
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

	//sphrayderv
	//lzw
	
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
	get_vars(fp_spc, "kmax ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &kmax);

	//runlsqr
	get_vars(fp_spc, "damper", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%f", &damper);
	}
	get_vars(fp_spc, "intlim", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%d", &intlims);
	}

	//runlsqr
	// c--Optionally read in some variables
	get_vars(fp_spc, "limitu ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &limitu);
	if (limitu != 0 && limitu != 1) {
		limitu = 0;
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

}
void read_files(char *spec_file){
	char *files[MUSTF] = { "oldvfil\0", "onedfil\0", //make1d 
						   "tgrdfil\0", "finevel\0", //c2f
						   "leqsfil\0", "fsumfil\0", "outlfil\0", "fhedfil\0", "fdatfil\0", // sphfdloc
						   "stafile\0", "locdfil\0", "telrerr\0", "dtdsfil\0", "resfile\0", "hitfile\0", "dtdhfil\0", "bookfil\0", "sclefil\0", //sphrayderv
						   "nmodfil\0", "fresfil\0", //runlsqr 
						   "fmodfil\0" //makenewmod
	};
	char *file_list[MUSTF] = { oldvfil, onedfil, tgrdfil, finevel, leqsfil, fsumfil, outlfil,
             				   fhedfil, fdatfil, stafile, locdfil, telrerr, dtdsfil, resfile, 
							   hitfile, dtdhfil, bookfil, sclefil, nmodfil, fresfil,  fmodfil };
	char pval[MAXSTRLEN + 1];
    int len, ierr;
	FILE *fp_spc;

    fp_spc = fopen(spec_file, "r");
	if (!fp_spc) {
		printf("(Error in read_spec.c)read fp_spc file error.\n");
		assert(0);
	}

	int lenf1, lenf2, i;
	for (i = 0; i < MUSTF; i++) {
		get_vars(fp_spc, files[i], pval, &len, &ierr);
		if (ierr == 1) {
			printf("Error trying to read filename %s", files[i]);
			assert(0);
		}
		sscanf(pval, "%s", file_list[i]);
	}
}

void read_error(char *name, char *type, FILE *fp_spc){
	printf("Error trying to read %s %s\n", type, name);
	fclose(fp_spc);

}