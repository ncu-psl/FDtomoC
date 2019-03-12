#ifndef EARTHQUAKE_FILE_DELIMITER
#define EARTHQUAKE_FILE_DELIMITER

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <string.h>
#include <sys/stat.h>
#ifdef _WIN32
	#include "windows/dirent.h"
#elif __linux__
	#include <dirent.h>
#endif

#include "parameter.h"
void clean_dir(char *);
void output_data(char *, int , char inp[maxobs][500], int);

int earthquake_file_delimiter(char *, char *);
#endif // EARTHQUAKE_FILE_DELIMITER
