#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <dirent.h>

#include "../include/c_parameter.h"
#include "../include/c_earthquake_file_delimiter.h"

void clean_dir(char *eqkdir) {
	struct dirent *de;
	DIR *dr = opendir(eqkdir);
	if (!dr) {
		printf("Could not open current directory: %s\n", eqkdir);
		return ;
	}
	while ((de = readdir(dr)) != NULL) {
		char filename[100];
		strcpy(filename,de->d_name);
		if(strcmp(filename,".")!=0 && strcmp(filename,"..")!=0) {
			strcpy(filename, eqkdir);
			strcat(filename,de->d_name);
			if(remove(filename)!=0) {
				printf("file can't remove: %s\n", filename);
				assert(0);
			}
		}
	}
	closedir(dr);
}

void output_data(char *eqkdir, int file_id, char inp[maxobs][500], int rows) {
	char filename[100];
	strcpy(filename, eqkdir);
	char tmp[100];
	sprintf(tmp, "%d%s", file_id, ".eqk");
	strcat(filename, tmp);
	FILE *fp = fopen(filename, "w");
	if(!fp) {
		printf("%s\n", filename);
		assert(0);
	}
	for (int i = 0; i < rows; i++) {
		fputs(inp[i], fp);
	}
	fclose(fp);
}

int earthquake_file_delimiter(char *filename, char *eqkdir) {
	mkdir(eqkdir, 0755);
	clean_dir(eqkdir);

	FILE *fp = fopen(filename, "r");
	if(!fp) {
		printf("open file error\n");
		printf("filename=%s\n", filename);
		assert(0);
	}
	int file_id = 1;
	for (file_id = 1; !feof(fp); file_id++) {
		char inp[maxsta][500];
		memset(inp, 0, sizeof(inp));
		for (int row = 0; row < maxobs; row++) {
			if (fgets(inp[row], 500, fp) == NULL)
				break;
			if (inp[row][0] == '\n') {
				output_data(eqkdir, file_id, inp, row);
				break;
			}
		}
	}
	fclose(fp);
	return file_id - 2;
}
