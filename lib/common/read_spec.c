#include "common/read_spec.h"


#define MUSTV  4
#define MUSTF  22
#define MAXSTRLEN 132

void read_variables(char *spec_file, SPEC *spec){
    char *mvals[MUSTV] = { "nxc\0", "nyc\0", "nzc\0", "h\0" };
	char pval[MAXSTRLEN + 1];
    int len, ierr;
	sscanf(spec_file, "%s", spec->spec_file);
    FILE *fp_spc;

    fp_spc = fopen(spec->spec_file, "r");
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
			sscanf(pval, "%d", &spec->nxc);
		} else if (i == 1) {
			sscanf(pval, "%d", &spec->nyc);
		} else if (i == 2) {
			sscanf(pval, "%d", &spec->nzc);
		} else if (i == 3) {
			sscanf(pval, "%lf", &spec->h);
		}
	}

//----dimension check
	if (spec->nxc > nxcm) {
		printf("nxc is too large, maximum is: %d\n", nxcm);
		assert(spec->nxc <= nxcm);
	} else if (spec->nyc > nycm) {
		printf("nyc is too large, maximum is: %d\n", nycm);
		assert(spec->nyc <= nycm);
	} else if (spec->nzc > nzcm) {
		printf("nzc is too large, maximum is: %d\n", nzcm);
		assert(spec->nzc <= nzcm);
	}

	//--Optionally read in some variables
//---Coordinate origin (used in header)
	get_vars(fp_spc, "x0 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &spec->x0);
	get_vars(fp_spc, "y0 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &spec->y[0]);
	get_vars(fp_spc, "z0 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &spec->z0);
	get_vars(fp_spc, "clat ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &spec->clat);
	get_vars(fp_spc, "clon ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &spec->clon);
	get_vars(fp_spc, "cz ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &spec->cz);
	get_vars(fp_spc, "azmod ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &spec->az);
	get_vars(fp_spc, "azmod ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &spec->azmod);
	get_vars(fp_spc, "ittnum ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->ittnum);
	get_vars(fp_spc, "df ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &spec->df);
	get_vars(fp_spc, "dq ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &spec->dq);

//----flatness, Vs, and sph  flags
	get_vars(fp_spc, "flat ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->iflat);

	get_vars(fp_spc, "vs1d ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->vs1d);

	get_vars(fp_spc, "sph ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->isph);

	//sphfdloc
	get_vars(fp_spc, "iread ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%d", &spec->iread);
	}
	if (spec->iread != 0 && spec->iread != 1) {
		spec->iread = 0;
	}
	get_vars(fp_spc, "ivs ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%d", &spec->ivs);
	}
	if (spec->ivs != 0 && spec->ivs != 1) {
		spec->ivs = 0;
	}
	get_vars(fp_spc, "vpvs ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &spec->vpvs);
	get_vars(fp_spc, "ivpvs ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%d", &spec->ivpvs);
	}
	if (spec->ivpvs != 0 && spec->ivpvs != 1) {
		spec->ivpvs = 0;
	}

	get_vars(fp_spc, "timedir ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", spec->timedir);
	get_vars(fp_spc, "eqkdir ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%s", spec->eqkdir);
	if(spec->eqkdir[strlen(spec->eqkdir) - 1] != '/'){
		spec->eqkdir[strlen(spec->eqkdir)] = '/';
		spec->eqkdir[strlen(spec->eqkdir) + 1] = '\0';
	}

	get_vars(fp_spc, "nthres ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->nthres);
	get_vars(fp_spc, "resthres ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &spec->resthres);
	get_vars(fp_spc, "resthrep ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &spec->resthrep);
	get_vars(fp_spc, "stdmax ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%lf", &spec->stdmax);
	get_vars(fp_spc, "kmin ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%d", &spec->kmin);
		spec->kmin--;
	}

//-----grid search control
	get_vars(fp_spc, "ndiv ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->ndiv);
	if (spec->ndiv <= 0)
		spec->ndiv = 1;
	get_vars(fp_spc, "ndiv2 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->ndiv2);
	if (spec->ndiv2 <= 0)
		spec->ndiv2 = 1;

	get_vars(fp_spc, "total_earthquakes ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->total_earthquakes);

	//sphrayderv
	//lzw
	
	get_vars(fp_spc, "vpvsscale ", pval, &len, &ierr);
	if (ierr == 0 && spec->ivpvs == 1) {
		sscanf(pval, "%f", &spec->vpvsscale);
	}

	get_vars(fp_spc, "dmean ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->idmean);
	if (spec->idmean != 0 && spec->idmean != 1) {
		spec->idmean = 0;
	}
	get_vars(fp_spc, "iray ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->iray);
	if (spec->iray != 0 && spec->iray != 1) {
		spec->iray = 0;
	}
	get_vars(fp_spc, "iraystat ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->iraystat);
	if (spec->iraystat != 0 && spec->iraystat != 1) {
		spec->iraystat = 0;
	}
	get_vars(fp_spc, "idatout ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->idatout);
	if (spec->idatout != 0 && spec->idatout != 1) {
		spec->idatout = 0;
	}
	get_vars(fp_spc, "nomat ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->nomat);
	if (spec->nomat != 0 && spec->nomat != 1) {
		spec->nomat = 0;
	}
	get_vars(fp_spc, "istacor ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->istacor);
	if (spec->istacor != 0 && spec->istacor != 1) {
		spec->istacor = 0;
	}
	get_vars(fp_spc, "doshot ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->idoshot);
	if (spec->idoshot != 0 && spec->idoshot != 1) {
		spec->idoshot = 0;
	}
	get_vars(fp_spc, "dotel ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->idotel);
	if (spec->idotel != 0 && spec->idotel != 1) {
		spec->idotel = 0;
	}
	get_vars(fp_spc, "resflag ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &spec->resflag);
	get_vars(fp_spc, "kmax ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->kmax);

	//runlsqr
	get_vars(fp_spc, "damper", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%f", &spec->damper);
	}
	get_vars(fp_spc, "intlim", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%d", &spec->intlims);
	}

	//runlsqr
	// c--Optionally read in some variables
	get_vars(fp_spc, "limitu ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->limitu);
	if (spec->limitu != 0 && spec->limitu != 1) {
		spec->limitu = 0;
	}
	get_vars(fp_spc, "mavx ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->mavx);
	get_vars(fp_spc, "mavy ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->mavy);
	get_vars(fp_spc, "mavz ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->mavz);
	get_vars(fp_spc, "nsmooth ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->nsmooth);
	get_vars(fp_spc, "dvperc ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &spec->dvperc);
	get_vars(fp_spc, "ipscflg ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->ipscflg);
	get_vars(fp_spc, "pertscl ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &spec->pertscl);
	get_vars(fp_spc, "do1d ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &spec->ido1d);

}
void read_files(char *spec_file, SPEC *file_identifier){
	char *files[MUSTF] = { "oldvfil\0", "onedfil\0", //make1d 
						   "tgrdfil\0", "finevel\0", //c2f
						   "leqsfil\0", "fsumfil\0", "outlfil\0", "fhedfil\0", "fdatfil\0", // sphfdloc
						   "stafile\0", "locdfil\0", "telrerr\0", "dtdsfil\0", "resfile\0", "hitfile\0", "dtdhfil\0", "bookfil\0", "sclefil\0", //sphrayderv
						   "nmodfil\0", "fresfil\0", //runlsqr 
						   "fmodfil\0", //makenewmod
						   "parlist\0" //sphfd
	};
	char *file_list[MUSTF] = { file_identifier->oldvfil, file_identifier->onedfil, file_identifier->tgrdfil, file_identifier->finevel, file_identifier->leqsfil, file_identifier->fsumfil, file_identifier->outlfil,
             				   file_identifier->fhedfil, file_identifier->fdatfil, file_identifier->stafile, file_identifier->locdfil, file_identifier->telrerr, file_identifier->dtdsfil, file_identifier->resfile, 
							   file_identifier->hitfile, file_identifier->dtdhfil, file_identifier->bookfil, file_identifier->sclefil, file_identifier->nmodfil, file_identifier->fresfil,  file_identifier->fmodfil,
							   file_identifier->parlist };
	char pval[MAXSTRLEN + 1];
    int len, ierr;
	FILE *fp_spc;

    fp_spc = fopen(spec_file, "r");
	if (!fp_spc) {
		printf("(Error in read_spec.c)read fp_spc file error.\n");
		assert(0);
	}

	for (int i = 0; i < MUSTF; i++) {
		get_vars(fp_spc, files[i], pval, &len, &ierr);
		if (ierr == 1) {
			printf("Error trying to read filename %s", files[i]);
			assert(0);
		}
		sscanf(pval, "%s", file_list[i]);
	}

	char *files_optional[MUSTF] = { "telefil\0",  "pbasfil\0", "sbasfil\0", "shotfil\0","elipfil\0", "raystat\0",  
									"dotfile\0", "headfil\0", "entfile\0", "stcfile\0" //sphrayderv
	};
	char *file_opt_list[MUSTF] = { file_identifier->telefil, file_identifier->pbasfil, file_identifier->sbasfil, file_identifier->shotfil, file_identifier->elipfil,
									file_identifier->raystat, file_identifier->dotfile, file_identifier->headfil, file_identifier->entfile, file_identifier->stcfile

	};

	for (int i = 0; i < 10; i++) {
		get_vars(fp_spc, files_optional[i], pval, &len, &ierr);
		if (ierr == 1) {
			printf("Error trying to read filename %s", files[i]);
			assert(0);
		}
		sscanf(pval, "%s", file_opt_list[i]);
	}

}

void read_grid(char *spec_file, SPEC *spec){	
	FILE *fp_spc;

    fp_spc = fopen(spec_file, "r");
	if (!fp_spc) {
		printf("(Error in read_spec.c)read fp_spc file error.\n");
		assert(0);
	}

	char aline[MAXSTRLEN + 1], varname[MAXSTRLEN + 1], parval[MAXSTRLEN + 1];
	int len, ierr;
	int ib = 0, ie = 0, lenv = 0, nvl = 0;
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
	sscanf(parval, "%d", &spec->igridx[0]);
	int k;
	for (k = 1; k < spec->nxc; k++) {
		ib = ie;
		get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
		sscanf(parval, "%d", &spec->igridx[k]);
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
	sscanf(parval, "%d", &spec->igridy[0]);
	for (k = 1; k < spec->nyc; k++) {
		ib = ie;
		get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
		sscanf(parval, "%d", &spec->igridy[k]);
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
		sscanf(parval, "%d", &spec->igridz[0]);
		for (k = 1; k < spec->nzc; k++) {
			ib = ie;
			get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
			sscanf(parval, "%d", &spec->igridz[k]);
		}
	}
	fclose(fp_spc);

}

void read_error(char *name, char *type, FILE *fp_spc){
	printf("Error trying to read %s %s\n", type, name);
	fclose(fp_spc);

}