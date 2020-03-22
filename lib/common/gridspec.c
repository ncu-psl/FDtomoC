#include <assert.h>
#include "common/geographic_method.h"
#include "common/read_spec.h"
#include "common/gridspec.h"
void moveGrid(GRID *grid){

    grid->nx = 1;
	grid->ny = 1;
	grid->nz = 1;
	grid->gx[0] = grid->x0;
	grid->gy[0] = grid->y[0];
	grid->gz[0] = grid->z0;

	for (int i = 1; i < grid->nxc; i++) {
		grid->nx = grid->nx + grid->igridx[i - 1];
		grid->gx[i] = grid->gx[i - 1] + grid->xSpace * grid->igridx[i - 1];
	}

	for (int i = 1; i < grid->nyc; i++) {
		grid->ny = grid->ny + grid->igridy[i - 1];
		grid->gy[i] = grid->gy[i - 1] + grid->ySpace * grid->igridy[i - 1];
	}

	for (int i = 1; i < grid->nzc; i++) {
		grid->nz = grid->nz + grid->igridz[i - 1];
		grid->gz[i] = grid->gz[i - 1] + grid->h * grid->igridz[i - 1];
	}
}


void coordinateCheck(GRID *grid, int isph, int isElevation){
    double z0r;
    double rearth = 6371.0f, degrad = 0.017453292f, hpi = 1.570796f;

    if (isph == 1) {
		grid->y[0] = grid->y00 * degrad;
//		y00 = hpi - glat(y00)
		if(isElevation){
			grid->y[0] = hpi - glath(grid->y[0], grid->z0, &z0r);
	} else {
			grid->y[0] = hpi - glat(grid->y[0]);
		}
		grid->x0 = grid->x00 * degrad; 
		grid->dq = grid->h / rearth;
		grid->df = fabs(grid->h / (rearth * sin(grid->y[0])));
		grid->dy = grid->dq / degrad;
		grid->dx = grid->df / degrad;
		grid->xSpace = grid->df;
		grid->ySpace = grid->dq;
	} else {
		grid->y[0] = grid->y00;
		grid->x0 = grid->x00;
		grid->dx = grid->h;
		grid->dy = grid->h;
		grid->xSpace = grid->h;
		grid->ySpace = grid->h;
	}

}
void dimensionCheck(GRID grid){

	if (grid.nx > nxm) {
		printf(" nx is too large, maximum is: %d", nxm);
		assert(!(grid.nx > nxm));
	}
	if (grid.ny > nym) {
		printf(" ny is too large, maximum is: %d", nym);
		assert(!(grid.ny > nym));
	}
	if (grid.nz > nzm) {
		printf(" nz is too large, maximum is: %d", nzm);
		assert(!(grid.nz > nzm));
	}
}