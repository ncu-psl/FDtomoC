#ifndef STRING_PROCESS
#define STRING_PROCESS

#include <string.h>
#include <stdio.h>
#include <assert.h>

void d_blank(char *, int *);
char *trim(char *);
char *strapp(char *, int *, const char *);
char *dtoa(char *, double, int);
void hdr_appender(char *, int, const char *, const char *, const char *, const char *, const char *, const char *);

#endif // STRING_PROCESS