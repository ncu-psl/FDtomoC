/*
c---General code to generate a coarse mesh 1D model from an
c	ASCII file specification and a spec file.
c
c	Author:	Steve Roecker, RPI
c
c	Version 2004.0909
c
c	Input files:
c
c	 Name   Uname   U#  Comments
c       onedfil lunone  14  ASCII file describing 1D model
c
c	Output files:
c
c	 Name   Uname   U#  Comments
c       oldvfil luncor  12  Coarse grid wavespeed model file
c
c	oldvfil is specified in the Spec file.
c
c	The format of the  1D input file should be:
c
c	<Descriptive Header>
c	Z1	Vp1	Vs1	C/I
c	Z2	Vp2	Vs2	C/I
c	...
c
c	where
c	Vp	is the P wavespeed
c	Vs	is the S wavespeed (or Vp/Vs if vs1d = 0)
c	Z	is the starting (upper) depth
c	C/I	is the interpolation specifier.  Set
c		to C for constant speed to next depth, and to
c		I to interpolate to the next depth.  If
c		the last wavespeed in the list is set to I, then
c		speeds at deeper depths are extrapolated using the
c		last available gradient (the difference between the
c		last two speeds divided by the difference between the
c		last two depths).
c
c	Variables required by this program:
c
c	nxc	Number of coarse grid points in the x direction
c	nyc	Number of coarse grid points in the y direction
c	nzc	Number of coarse grid points in the z direction
c	h	Find grid spacing
c
c	Optional variables that can be usefully set by the user:
c
c	igridx	Array of spaces between fine and coarse grid in x direction
c	igridy	Array of spaces between fine and coarse grid in y direction
c	igridz	Array of spaces between fine and coarse grid in z direction
c	flat	= 1 to flatten depths and speeds, 0 to leave as is. NB:  Because the tomo
c			sofware uses a regular grid, we presume that the input grid description 
c			is already flattened.  This program starts by "unflattening" these
c			depths prior to looking up the appropriate wavespeed from the table.
c			The output is flattened wavespeeds.
c	vs1d	= 1 to interpret Vs as shear wave speed; = 0 to interpret as Vp/Vs
c	isph	= 1 to use a spherical coordinate system; =0 for cartesian.
c		    note that the only practical difference is that the definitions of dx and dy
c		    change so that dx = df and dy = df.
*/

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "common/environment_setting.h"
#include "common/parseprogs.h"
#include "common/parameter.h"
#include "common/gridspec.h"
#include "common/geographic_method.h"
#include "common/string_process.h"
#include "common/shared_variables.h"
#include "common/vhead.h"
#include "common/interpolation.h"
#include "FDtomo/make1d.h"
#define MAXSTRLEN 132
//---number of 4 byte words in the header
#define nhbyte 58 * 4

float flatvel(float, float);
float uflatz(float);
float flatz(float);
char * dtoa(char *, double, int);
velocity3D create3DModel(SPEC spec, velocity1D model) {
	float *gz = spec.grid.gz;
	//---unflatten the depths if required
	if (spec.iflat == 1) {
		for (int i = 0; i < spec.grid.nzc; i++) {
			gz[i] = uflatz(gz[i]);
		}
	}
//----generate the model
	float *vp = (float *)malloc(sizeof(float) * spec.grid.nzc);
	float *vs = (float *)malloc(sizeof(float) * spec.grid.nzc);
	vp = linear_interpolation_array(spec.grid.gz, model.z, model.vp, model.nl, spec.grid.nzc, model.terp);
	vs = linear_interpolation_array(spec.grid.gz, model.z, model.vs, model.nl, spec.grid.nzc, model.terp);
	velocity3D model3D = generate3DModel(vp, vs, spec.grid);
	return model3D;
}

float flatvel(float v, float z) {
// --- does earth-flattening correction for velocity, at depth z.
//     assumes z is the flat-earth depth (already corrected). 12/87 gaa.
//     a different transform is done for blocks, which have integrated
//     average velocities
	float r = 6371.00f;
	return expf(z / r) * v;
}

float flatz(float z) {
// --- does earth-flattening correction from depth in spherical earth to
//     depth in flat-earth.  this preserves travel-times if the velocity
//     correction is also used.  see chapman (1973).  gaa 12/87.
	float r = 6371.00f;
	return r * logf(r / (r - z)) / logf(expf(1));
}

float uflatz(float z) {
// --- undoes earth-flattening correction from depth in spherical earth to
//     depth in flat-earth.

	float r = 6371.00f;
	return r * (1. - expf(-z / r));
}

int OUTPUT_MAKE1D(MAKE1D_DATA *maked1d_data, SPEC spec){
	FILE *fp_cor;
	fp_cor = fopen(spec.oldvfil, "wb");
	if (!fp_cor) {
		printf("file create error: %s\n", spec.oldvfil);
		assert(0);
	}
	
	int nxc = spec.grid.nxc, nyc = spec.grid.nyc, nzc = spec.grid.nzc;
	int nxyc = nxc * nyc;
	int nxyzc = nxyc * nzc;
	int nxyzc2 = nxyzc * 2;
	fwrite(&maked1d_data->head, sizeof(char), nhbyte, fp_cor);
	fwrite(maked1d_data->igridx, sizeof(maked1d_data->igridx[0]), nxc - 1, fp_cor);
	fwrite(maked1d_data->igridy, sizeof(maked1d_data->igridy[0]), nyc - 1, fp_cor);
	fwrite(maked1d_data->igridz, sizeof(maked1d_data->igridz[0]), nzc - 1, fp_cor);
	fwrite(maked1d_data->vsave, sizeof(maked1d_data->vsave[0]), nxyzc2, fp_cor);
	return 0;
}

int LOG_MAKE1D(SPEC spec){
	char VERSION[10] = "2004.0909\0";
	int nxyc = spec.grid.nxc * spec.grid.nyc;
	int nxyzc = nxyc * spec.grid.nzc;
	int nxyzc2 = nxyzc * 2;

	int nxy = spec.grid.nx * spec.grid.ny;
	int nxyz = nxy * spec.grid.nz;
	
	int lenhead = nhbyte;
	int lengrd = 4 * (spec.grid.nxc + spec.grid.nyc + spec.grid.nzc - 3);
	int lenrec = lenhead + lengrd + 4 * nxyzc2;

	FILE *fp_log;

	fp_log = fopen("make1d.log", "w");
	if(!fp_log) {
		printf("(Error in make1d.c)create fp_log file error.\n");
		assert(fp_log);
	}

	fprintf(fp_log, "  \n");
	fprintf(fp_log,
			" *************************************************************** \n");

	fprintf(fp_log, "          Parameters Set For This Run of make1d.f\n");

	fprintf(fp_log, "  \n");
	fprintf(fp_log, " VERSION: %s\n", VERSION);
	fprintf(fp_log, "  \n");
	fprintf(fp_log, " Current parameter specification file: %-40s\n",
			spec.spec_file);
	fprintf(fp_log, "  \n");
	char out_str[MAXSTRLEN];
	dtoa(out_str, spec.clat, 18);
	fprintf(fp_log, "  Latitude origin  (clat):   %s     \n", out_str);
	dtoa(out_str, spec.clon, 18);
	fprintf(fp_log, "  Longitude origin (clon):   %s     \n", out_str);
	dtoa(out_str, spec.cz, 18);
	fprintf(fp_log, "  Depth of  origin (cz)  :   %s     \n", out_str);
	dtoa(out_str, spec.az, 10);
	fprintf(fp_log, "  Clockwise rotation (az):   %s    \n", out_str);
	fprintf(fp_log, "  \n");
	dtoa(out_str, spec.grid.x0, 18);
	fprintf(fp_log, "  Cartesian X origin (x0):   %s     \n", out_str);
	dtoa(out_str, spec.grid.y[0], 18);
	fprintf(fp_log, "  Cartesian Y origin (y0):   %s     \n", out_str);
	dtoa(out_str, spec.grid.z0, 18);
	fprintf(fp_log, "  Cartesian Z origin (z0):   %s     \n", out_str);
	if (spec.isph == 0) {
		fprintf(fp_log, "  Coordinate system is CARTESIAN \n");
		fprintf(fp_log, "  Fine grid spacing: %lf\n", spec.grid.h);
	} else {
		fprintf(fp_log, "  Coordinate system is SPHERICAL \n");
		fprintf(fp_log, "  Fine Longitude spacing (df):   %.17E\n", spec.grid.dx);
		fprintf(fp_log, "  Fine Latidtude spacing (dq):   %.17E\n", spec.grid.dy);
		dtoa(out_str, spec.grid.h, 18);
		fprintf(fp_log, "  Fine Radial spacing    (dz):   %s     \n", out_str);
	}
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of X coarse grid nodes: %12d\n", spec.grid.nxc);
	fprintf(fp_log, "  X coarse grid node spacing: \n");
	for (int iii = 0; iii < spec.grid.nxc - 1; iii++) {
		fprintf(fp_log, "% 4d", spec.grid.igridx[iii]);
		if (iii % 10 == 9) {
			fprintf(fp_log, "\n");
		}
	}
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of Y coarse grid nodes: %12d\n", spec.grid.nyc);
	fprintf(fp_log, "  Y coarse grid node spacing: \n");
	for (int iii = 0; iii < spec.grid.nyc - 1; iii++) {
		fprintf(fp_log, "% 4d", spec.grid.igridy[iii]);
		if (iii % 10 == 9) {
			fprintf(fp_log, "\n");
		}
	}
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of Z coarse grid nodes: %12d\n", spec.grid.nzc);
	fprintf(fp_log, "  Z coarse grid node spacing: \n");
	for (int iii = 0; iii < spec.grid.nzc - 1; iii++) {
		fprintf(fp_log, "% 4d", spec.grid.igridz[iii]);
		if (iii % 10 == 9) {
			fprintf(fp_log, "\n");
		}
	}
	fprintf(fp_log, "\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Number of X fine grid nodes: %12d\n", spec.grid.nx);
	fprintf(fp_log, "  Number of Y fine grid nodes: %12d\n", spec.grid.ny);
	fprintf(fp_log, "  Number of Z fine grid nodes: %12d\n", spec.grid.nz);
	fprintf(fp_log, " \n");
	fprintf(fp_log, "  Total Number of coarse grid nodes:%13d\n", nxyzc);
	fprintf(fp_log, "  Total Number of fine grid nodes:%13d\n", nxyz);
	fprintf(fp_log, " \n");
	if (spec.iflat == 1) {
		fprintf(fp_log, "  Speeds and Depths are flattened (iflat = 1) \n");
	} else {
		fprintf(fp_log, "  Speeds and Depths are not flattened (iflat = 0) \n");
	}
	if (spec.vs1d == 1) {
		fprintf(fp_log, "  Vs column treated as S wavespeed (vs1d = 1) \n");
	} else {
		fprintf(fp_log, "  Vs column treated as Vp/Vs (vs1d = 0) \n");
	}

	fprintf(fp_log, " \n");

	fprintf(fp_log, "  Length of header  (bytes): %12d\n", lenhead + lengrd);
	fprintf(fp_log, "  Length of outfile (bytes): %12d\n", lenrec);

	fprintf(fp_log, " \n");
	fprintf(fp_log, " Input file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " One-D ASCII model file: %-40s\n", spec.onedfil);
	fprintf(fp_log, " \n");
	fprintf(fp_log, " Output file attachments:\n");
	fprintf(fp_log, " \n");
	fprintf(fp_log, " New Wavespeed model file: %-40s\n", spec.oldvfil);
	fprintf(fp_log, " \n");

	fclose(fp_log);
	return 0;

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
