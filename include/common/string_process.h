#ifndef STRING_PROCESS
#define STRING_PROCESS

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

void d_blank(char *, int *);
char *trim(char *);
size_t trimwhitespace(char *, size_t, const char *);
char *strapp(char *, int *, const char *);
char *dtoa(char *, double, int);
void hdr_appender(char *, int, const char *, const char *, const char *, const char *, const char *, const char *);

#endif // STRING_PROCESS