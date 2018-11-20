#include "../include/c_string_process.h"

void d_blank(char *line, int *len) {
	int len_line = strlen(line);
	char tmp[len_line];
	int j = 0;
	for (int i = 0; i < len_line; i++)
		if (line[i] != ' ') {
			tmp[j] = line[i];
			j++;
		}
	tmp[j] = '\0';
	int len_tmp = strlen(tmp);
	for (int i = 0; i < len_tmp; i++)
		line[i] = tmp[i];
	for (int i = len_tmp; i < len_line; i++)
		line[i] = '\0';
	*len = len_tmp;
}

void trim(char *str) {
	char *end;

	// Trim leading space
	while (((unsigned char) *str) == ' ')
		str++;

	if (*str == 0)  // All spaces?
		return;  // str;

	// Trim trailing space
	end = str + strlen(str) - 1;
	while (end > str && ((unsigned char) *end) == ' ')
		end--;

	// Write new null terminator
	*(end + 1) = 0;

	//return str;
}

char *strapp(char *dest, int *end, const char *src)
{
	if ( *end >= 0 && dest && src )
	{
		char *p = dest + *end;
		while ((*p++ = *src++ ))
			(*end)++;
	}
	return dest;
}

char * dtoa(char * str_output, double inp, int len) {
	double num = inp;
	if (inp < 0)
		num = -inp;

	snprintf(str_output, 129, " %.131lf", num);
	str_output[0] = ' ';
	if (inp < 0)
		str_output[0] = '-';
	str_output[len + 1] = '\0';

	//if(debug_print)
	//printf("%s\n",str_output);
	return str_output;
	/*
	 char *ptr = strchr(a, '.');
	 if(ptr) {
	 int index = ptr - a;
	 }else {

	 }
	 */
}

void hdr_appender(char *hdr, int nhbyte, const char *head, const char *type, const char *syst, const char *quant, const char *flatten, const char *hcomm) {
	hdr[0] = '\0';
	sprintf(hdr, "%4s%4s%4s%4s%4s", head, type, syst, quant, flatten);
	if(strlen(hcomm) > nhbyte - strlen(hdr)) {
		printf("string length of hcomm is too long: strlen(hcomm) = %lu \n", strlen(hcomm));
		assert(0);
	}
	strcat(hdr, hcomm);
	for (int i = strlen(hdr); i < nhbyte; i++)
		hdr[i] = ' ';
	hdr[nhbyte] = '\0';
}
