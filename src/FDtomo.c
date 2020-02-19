#include "FDtomo/make1d.h"
#include "FDtomo/c2f.h"
#include "FDtomo/sphfd.h"
#include "FDtomo/sphfdloc.h"
#include "FDtomo/sphrayderv.h"
#include "FDtomo/runlsqr.h"
#include "FDtomo/makenewmod.h"

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
	strcpy(spec.timedir, "./ttimes\0");
	strcpy(spec.eqkdir, "./eqkdir/\0");
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

	//sphrayderv initialization


	read_variables(argv[1], &spec);
	read_files(argv[1], &spec);
	read_grid(argv[1], &spec);
	MAKE1D_DATA *MAKE1D = make1d(spec);
	OUTPUT_MAKE1D(MAKE1D, spec);
	C2F_DATA *C2F = c2f(spec, MAKE1D);
	OUTPUT_C2F(C2F, spec);
	SPHFD_DATA **SPHFD = sphfd(argc, argv, spec, C2F);
	for (int i = 0; i < num_parfiles; i++)
		OUTPUT_SPHFD(SPHFD[i], spec);
	SPHFDLOC_DATA **SPHFDLOC = sphfdloc(spec, SPHFD);
	SPHRAYDERV_DATA *SPHRAYDERV = sphrayderv(spec, SPHFDLOC);
	RUNLSQR_DATA *RUNLSQR = runlsqr(spec, SPHRAYDERV);
	makenewmod(spec, RUNLSQR);
	return 0;

}
