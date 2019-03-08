#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int get_sta_chan(char *, char *, char *, int *);
void get_vars(FILE*, char *, char *, int*, int*);
void get_field(FILE *, char *, int, int *, char *, int *, int *);
void get_line(FILE*, char*, int *);
int get_var_somp(int, char*, char*);

int main(int argc, char** args) {
	int nx = 301, ny = 301, nz = 61, reverse = 50;
	char velfine[10]="VP.mod";
	FILE *fp_spc = fopen(args[1], "r");
	if(!fp_spc) {
		printf("Can not open file: %s\n", args[1]);
		assert(0);
	}
	FILE *fp_sta = fopen(args[2], "r");
	if(!fp_sta) {
		printf("Can not open file: %s\n", args[2]);
		assert(0);
	}

	int len = 0, ierr = 0;
	char pval[100];

	float x0 = 0, y0 = 0, z0 = 0, h = 0;
	get_vars(fp_spc, "x0 ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%f", &x0);
	}
	get_vars(fp_spc, "y0 ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%f", &y0);
	}
	get_vars(fp_spc, "z0 ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%f", &z0);
	}
	get_vars(fp_spc, "h ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%f", &h);
	}

	int sph = 0;
	get_vars(fp_spc, "sph ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%f", &sph);
	}

	char timedir[100];
	get_vars(fp_spc, "timedir ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%s", &timedir);
	}
	fclose(fp_spc);

	char sphlist[100] = "runsphfd01.sh";
	FILE* fp_list = fopen(sphlist, "w");
	if(!fp_list) {
		printf("Can not create file: runsphfd01\n");
		assert(0);
	}

	char str_inp[100];
	while(fgets(str_inp, sizeof(str_inp), fp_sta) != NULL) {
		if(str_inp[0] == '\n') {
			break;
		}

		int zs = 0;
		float stx = 0, sty = 0, fxs = 0, fys = 0, fzs = 0;
		char sta[100];
		if(sph) {
			sscanf(str_inp, "%f %f %d %s %f %f", &sty, &stx, &zs, sta, &fys, &fxs);
		} else {
			sscanf(str_inp, "%f %f %f %s", &fys, &fxs, &fzs, sta);
		}
		fzs = zs / 1000.;

		for(int j=0;j<2;j++) {
			char filename[100];
			strcpy(filename, sta);
			if(j==0) {
				strcat(filename,".par");
				fprintf(fp_list, "./sphfd par=../parfiles_P/%s.par\n", sta);
			} else {
				strcat(filename,".spar");
				fprintf(fp_list, "./sphfd par=../parfiles_P/%s.spar\n", sta);
			}

			FILE* fp_time = fopen(filename, "w");
			fprintf(fp_time, "fxs=%f\nfys=%f\nfzs=%f\n", fxs, fys, fzs);
			fprintf(fp_time, "nx=%f\nny=%f\nnz=%f\n", (float)nx, (float)ny, (float)nz);
			fprintf(fp_time, "x0=%f\ny0=%f\nz0=%f\nh=%f\n", x0, y0, z0, h);
			if(j == 0) {
				fprintf(fp_time, "timefile=%s.ptimes\n", sta);
				fprintf(fp_time, "velfile=VP.mod\n");
			} else {
				fprintf(fp_time, "timefile=%s.stimes\n", sta);
				fprintf(fp_time, "velfile=VS.mod\n");
			}
			fprintf(fp_time, "reverse=50\n");

			fclose(fp_time);
		}
	}
	fclose(fp_list);
	fclose(fp_sta);
	return 0;
}


#define MAXSTRLEN 132
/*
 c-- - recover station and channel name from PASSCAL file name
 */

void get_vars(FILE *fp, char *vname, char *parval, int *len, int *error) {
	if (strlen(vname) >= MAXSTRLEN) {
		printf("(Error in get_vars) vname's length is larger than 132: %s\n",
			vname);
		assert(0);
	}
	memset(parval, 0, strlen(parval) * sizeof(char));
	char aline[MAXSTRLEN];
	char varname[MAXSTRLEN];
	int ierr = 0, i = 0, lenp = (int)strlen(vname);
	*error = 1;

	//c---recover the length of vname
	while (i < (int)strlen(vname)) {
		if (vname[i] == ' ') {
			//go to a2
			lenp = i;
			break;
		}
		i++;
	}

	//a2: lenp = i;

	rewind(fp);

	int nline = 0;
a1: memset(aline, 0, sizeof(aline));
	get_line(fp, aline, &ierr);
	if (ierr == 1)
		goto a60;
	nline = nline + 1;
	if (ierr != 0)
		goto a1;

	//c------recover the variable name
		//int ib = 1, ie = 0, lenv = 0;
	int ib = 0, ie = 0, lenv = 0;
	memset(varname, 0, sizeof(varname));
	get_field(fp, aline, ib, &ie, varname, &lenv, &ierr);
	//printf("60 aline=%s varname=%s \n",aline,varname);
	if (ierr != 0)
		goto a15;
	//c-----see if this parameter name is valid
	if (lenv != lenp)
		goto a1;
	if (strncmp(varname, vname, lenv) != 0)
		goto a1;

	//c-----recover the variable setting
	ib = ie;
	/////////////////////
	int nvl;
	get_field(fp, aline, ib, &ie, parval, &nvl, &ierr);
	if (ierr != 0)
		goto a15;

	//c	write(*,*) ' j = ', j
	//c	write(*,*) ' parval:  ', parval(1:nvl)
	*len = nvl;
	*error = 0;
	return;

a15: printf(" get_vars Warning:  get_field error at line %d\n", nline);
	goto a1;

a60: *error = 1;
	//c	write(*,*) ' getvars warning:  EOF reached in pararmeter file'

	return;
}

void get_field(FILE *fp, char *aline, int ib, int *ie, char *field, int *len,
	int *ierr) {
	if (strlen(aline) >= MAXSTRLEN) {
		printf("(Error in get_field) aline's length is larger than %d: %s\n",
			MAXSTRLEN, aline);
		*ierr = 1;
		assert(0);
	}
	if (strlen(aline) == 0) {
		printf("(Error in get_field) aline's length is zero: %s\n", aline);
		*ierr = 1;
		assert(0);
	}

	char tab = '\t';
	int i, i1, nch;
a1: nch = 0;
	for (i = ib; i < MAXSTRLEN; i++) {
		//c----ignore leading tabs or spaces
		if (nch == 0 && (aline[i] == ' ' || aline[i] == tab))
			continue;
		//goto a10;
		if (aline[i] == ' ' || aline[i] == tab || aline[i] == '\0')
			goto a12;
		//c----test for a quoted string
		if (nch == 0 && aline[i] == '"')
			goto a14;
		//c----test for a continuation
		if (nch == 0 && aline[i] == '\\')
			goto a22;
		field[nch] = aline[i];
		nch = nch + 1;
		*len = nch;
	}
	printf(" Error:  EOL reached in get_field! \n");
	field[0] = '\0';
	*ierr = 1;
	return;

a12: *ierr = 0;
	*ie = i;
	field[nch] = '\0';
	return;

	//c----quoted string section
a14: i1 = i + 1;
	nch = 0;
	for (i = i1; i < MAXSTRLEN; i++) {
		if (aline[i] == '"')
			goto a12;
		field[nch] = aline[i];
		nch = nch + 1;
		*len = nch;
	}
	field[MAXSTRLEN - 1] = '\0';
	printf(" Error:  End of quoted string not found in getfield! \n");
	*ierr = 1;
	return;

	//c----continuation section
a22: get_line(fp, aline, ierr);
	if (*ierr == 1)
		goto a60;
	if (*ierr != 0)
		goto a50;
	ib = 1;
	goto a1;

a50: printf(" Read error in get_field! \n");
	field[0] = '\0';
	*ierr = 1;
	return;

a60: printf(" Error:  EOF encountered in get_field! \n");
	field[0] = '\0';
	*ierr = 1;

	return;
}

void get_line(FILE *fp, char *aline, int *ierr) {
	if (strlen(aline) >= MAXSTRLEN) {
		printf("(Error in get_line) aline's length is larger than %d: %s\n",
			MAXSTRLEN, aline);
		*ierr = 1;
		return;
	}
	aline[MAXSTRLEN - 1] = '\0';
	char tab = '\t';

a1: if (fgets(aline, MAXSTRLEN, fp) == NULL) {
	goto a60;
}
	else {
	//replace \n by \0
	if (strlen(aline) > 0) {
		aline[strlen(aline) - 1] = '\0';
	}
	else {
		printf("(Error in get_line) is an empty line\n");
	}
}

	//c------skip over comments
	if (aline[0] == '#' || aline[0] == '\0')
		goto a1;
	int i = 0;
	while (i < MAXSTRLEN) {
		if (aline[i] != ' ' && aline[i] != tab)
			goto a2;
		i++;
	}
	goto a1;
a2: *ierr = 0;
	return;

	// a50: printf(" Read error in get_line! \n");
	*ierr = 2;
	return;

a60: *ierr = 1;

	return;
}