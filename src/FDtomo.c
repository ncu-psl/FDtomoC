#include "FDtomo/make1d.h"
#include "common/shared_variables.h"

#include <stdio.h>
int main(int argc, char *argv[]){
	SPEC spec;
	//set default value
	//c2f initialization 
	spec.vs1d = 1;
	spec.iflat = 0;
	spec.isph = 0;

	//sphfdloc initialization 
	sscanf("./ttimes\0", "%s", spec.timedir);
	sscanf("./eqkdir/\0", "%s", spec.eqkdir);
	spec.iread = 0;
	spec.ivs = 1;
	spec.vpvs = 1.78;
	spec.nthres = 8;
	spec.resthres = .5;
	spec.resthrep = 5.0;
	spec.stdmax = 15.0;
	spec.kmin = 2;
	spec.x0 = 0;
	spec.y[0] = 0;
	spec.z0 = 0;
	spec.dq = 0;
	spec.df = 0;
	spec.ndiv = 20;
	spec.ndiv2 = 20;
	spec.ittnum = 1;
	spec.total_earthquakes = 0;


	read_variables(argv[1], &spec);
	read_files(argv[1], &spec);
	read_grid(argv[1], &spec);
	make1d(spec);
	c2f(spec);
	sphfd(argc, argv, spec);

	sphfdloc(spec);
	//sphrayderv(argv[1]);
	//runlsqr(argv[1]);
	//makenewmod(argv[1]);
	return 0;

}
