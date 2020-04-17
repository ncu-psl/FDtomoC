#include "common/interpolation.h"
#include "common/velocity_model.h"

velocity1D read_velocity1D(SPEC spec){
	velocity1D model;
	
	FILE *fp_one;
	fp_one = fopen(spec.onedfil, "r");
	if(!fp_one) {
        printf("Error on opening fp_one(%s)\n", spec.onedfil);
        assert(0);
	}

	char aline[MAXSTRLEN + 1];
	char pval[MAXSTRLEN + 1];

	double h;
    int len, ierr;
	int ib = 0, ie = 0, lenv = 0, nvl = 0;
//---skip over header
	get_line(fp_one, aline, &ierr);
	int nl = 0;
	a21: get_line(fp_one, aline, &ierr);
	if (ierr == 1){
		model.nl = nl;
		fclose(fp_one);
		return model;
	}
	if (ierr != 0)
		goto a21;
	ib = 0;
	//******************************************
	//      fp_one   or  fp_spc?
	//******************************************
	get_field(fp_one, aline, ib, &ie, pval, &nvl, &ierr);
	//read(pval(1:nvl),*, err=50) h
	sscanf(pval, "%lf", &h);
	ib = ie;
	get_field(fp_one, aline, ib, &ie, pval, &nvl, &ierr);
	float p = 0;
	sscanf(pval, "%f", &p);
	ib = ie;
	get_field(fp_one, aline, ib, &ie, pval, &nvl, &ierr);
	float s = 0;
	sscanf(pval, "%f", &s);
	ib = ie;
	get_field(fp_one, aline, ib, &ie, pval, &nvl, &ierr);

	if (nl >= MAX1D) {
		printf(" Too many layers; maximum now is %d", MAX1D);
		assert(!(nl >= MAX1D));
	}
	model.vp[nl] = p;
	if (spec.vs1d == 1) {
		model.vs[nl] = s;
	} else {
		model.vs[nl] = p / s;
	}
//****ONE-TIME CLUDGE TO FORCE A VP/VS
//   vp(nl,2) = p/1.78
//***********************************
	model.z[nl] = h;
	model.terp[nl] = pval[0];

	nl = nl + 1;
	goto a21;
	
}
velocity3D create3DModel(SPEC spec, velocity1D model) {
	float *gz = spec.grid.gz;
	//---unflatten the depths if required
    /*
	if (spec.iflat == 1) {
		for (int i = 0; i < spec.grid.nzc; i++) {
			gz[i] = uflatz(gz[i]);
		}
	}*/
//----generate the model
	float *vp = (float *)malloc(sizeof(float) * spec.grid.nzc);
	float *vs = (float *)malloc(sizeof(float) * spec.grid.nzc);
	vp = linear_interpolation_array(spec.grid.gz, model.z, model.vp, model.nl, spec.grid.nzc, model.terp);
	vs = linear_interpolation_array(spec.grid.gz, model.z, model.vs, model.nl, spec.grid.nzc, model.terp);
	velocity3D model3D = generate3DModel(vp, vs, spec.grid);
	return model3D;
}


velocity3D generate3DModel(float *vp, float *vs, GRID grid){
	velocity3D model3D;
	memcpy(&model3D.grid, &grid, sizeof(grid));
	int sizeOfGrid = grid.nxc * grid.nyc * grid.nzc;
	model3D.vp = (float *)malloc(sizeof(float) * sizeOfGrid);
	model3D.vs = (float *)malloc(sizeof(float) * sizeOfGrid);
	for (int i = 0; i < grid.nzc; i++){
		int nxyc = grid.nxc * grid.nyc;
		int ioff = nxyc * i;
		for(int j = 0; j < nxyc; j++){
			model3D.vp[ioff] = vp[i];
			model3D.vs[ioff] = vs[i];
		}
	}
	return model3D;
}
