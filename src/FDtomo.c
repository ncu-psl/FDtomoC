#include "FDtomo/make1d.h"
#include "common/read_spec.h"
#include <stdio.h>
int main(int argc, char *argv[]){
	read_variables(argv[1]);
	read_files(argv[1]);
	read_grid(argv[1]);
	make1d(argv[1]);
	c2f(argv[1]);
	sphfd(argc, argv, argv[1]);

	//sphfdloc(argv[1]);
	//sphrayderv(argv[1]);
	//runlsqr(argv[1]);
	//makenewmod(argv[1]);
	return 0;

}
