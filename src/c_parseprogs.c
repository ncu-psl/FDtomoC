#include "../include/c_parseprogs.h"
#include "../include/environment_setting.h"

#define MAXSTRLEN 132
/*

 c-- - recover station and channel name from PASSCAL file name

 */
//int get_sta_chan(const char *, char *, char *, int *);
//void get_vars(FILE*, char[132], char *, int*, int*);void get_field(FILE *, char *, int, int *, char *, int *, int *);
//void get_line(FILE*, char*, int *);
//int get_var_somp(int, char*, char*);

//not yet implemented
int get_sta_chan(char *_file_name, char *sta, char *chan, int *error) {
	return 0;
}

void get_vars(FILE *fp, char *vname, char *parval, int *len, int *error) {
	if (strlen(vname) >= MAXSTRLEN) {
		printf("(Error in get_vars) vname's length is larger than 132: %s\n",
				vname);
		assert(0);
	}
	memset(parval, 0, strlen(parval) * sizeof(char));
	char aline[MAXSTRLEN];
	char varname[MAXSTRLEN];
	int ierr = 0, i = 0, lenp = strlen(vname);
	*error = 1;

//c---recover the length of vname
	while (i < strlen(vname)) {
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
	field[MAXSTRLEN] = '\0';
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
		printf("(Error in get_field) aline's length is larger than %d: %s\n",
				MAXSTRLEN, aline);
		*ierr = 1;
		return;
	}
	aline[MAXSTRLEN] = '\0';
	char tab = '\t';

	a1: if (fgets(aline, MAXSTRLEN, fp) == NULL) {
		goto a60;
	} else {
		//replace \n by \0
		aline[strlen(aline) - 1] = '\0';
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

//not yet implemented
int get_var_somp(int lunspc, char *vname, char *parval) {
	return 0;
}
