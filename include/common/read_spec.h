#ifndef READ_SEPC
#define READ_SPEC
#include "parseprogs.h"
void read_variables(char *spec_file);
void read_files(char *spec_file);
void read_error(char *name, char *type, FILE *fp_spc);
#endif
