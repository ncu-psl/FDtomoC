#include "FDtomo/make1d.h"
#include "common/shared_variables.h"

#include <stdio.h>
int main(int argc, char *argv[]){
	SPEC spec;
	spec.vs1d = 1;
	spec.iflat = 0;
	spec.isph = 0;

	read_variables(argv[1], &spec);
	read_files(argv[1], &spec);
	read_grid(argv[1], &spec);
	make1d(spec);
	c2f(spec);
	//sphfd(argc, argv, argv[1]);

	//sphfdloc(argv[1]);
	//sphrayderv(argv[1]);
	//runlsqr(argv[1]);
	//makenewmod(argv[1]);
	return 0;

}
