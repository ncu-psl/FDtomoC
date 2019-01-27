#ifndef PARSEPROGS
#define PARSEPROGS

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <assert.h>

int get_sta_chan(char *, char *, char *, int *);
void get_vars(FILE*, char *, char *, int*, int*);
void get_field(FILE *, char *, int, int *, char *, int *, int *);
void get_line(FILE*, char*, int *);
int get_var_somp(int, char*, char*);

#endif // PARSEPROGS