#include "../include/sphrayderv/dfind.h"

void dfind(const double dd[3], double *d, int *md, double tolmin, double tolmax) {
	int mmin = 0, mmax = 0, mmid = 0;
	double dmin = 0, dmax = 0, dmid = 0;
	if (fabs(dd[0]) <= fabs(dd[1]) && fabs(dd[0]) <= fabs(dd[2])) {
		mmin = 0;
		if (fabs(dd[1]) <= fabs(dd[2])) {
			mmid = 1;
			mmax = 2;
		}
		else {
			mmid = 2;
			mmax = 1;
		}
	}
	else if (fabs(dd[1]) <= fabs(dd[2])) {
		mmin = 1;
		if (fabs(dd[0]) <= fabs(dd[2])) {
			mmid = 0;
			mmax = 2;
		}
		else {
			mmid = 2;
			mmax = 0;
		}
	}
	else {
		mmin = 2;
		if (fabs(dd[0]) <= fabs(dd[1])) {
			mmid = 0;
			mmax = 1;
		}
		else {
			mmid = 1;
			mmax = 0;
		}
	}

	dmin = dd[mmin];
	dmid = dd[mmid];
	dmax = dd[mmax];
	*d = dmin;
	*md = mmin;
	if (dmin < tolmin) {
		if (dmid < tolmin) {
			if (dmax < tolmax) {
				*d = dmax;
				*md = mmax;
			}
			else {
				*d = dmid;
				*md = mmid;
			}
		}
		else {
			if (dmid < tolmax) {
				*d = dmid;
				*md = mmid;
			}
		}
	}
}