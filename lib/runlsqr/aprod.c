// c
// c  aprod performs the following functions:
// c
// c      if mode = 1, set y = y + a*x
// c      if mode = 2, set x = x + a(transpose) * y
// c
// c  where a is a matrix stored by rows in the arrays a, ja, and na.
// c  in this example, a, ja, na are stored  in common.
// c

#include "runlsqr/aprod.h"

void aprod(int mode, int m, int n, float *x, float *y, int leniw, int lenrw, int *iw, float *rw, float *a, int *na,
		int *ja) {
	// x=(float *)malloc(sizeof(float)*n);
	// y=(float *)malloc(sizeof(float)*m);
	iw=(int *)malloc(sizeof(int)*leniw);
	rw=(float *)malloc(sizeof(float)*lenrw);

	if (mode == 1) {
		mode1(m, x, y, a, na, ja);
	} else if (mode == 2) {
		mode2(m, x, y, a, na, ja);
	} else {
		printf("mode = %d is unknow\n", mode);
		assert(0);
	}
	
	free(iw);
	free(rw);
}

// c--------------------------------------------------------------------------
// c mode = 1  -- set y = y + a*x.
// c--------------------------------------------------------------------------
void mode1(int m, float *x, float *y, float *a, int *na, int *ja) {
	int l1 = 0, l2 = 0;
	for (int i = 0; i < m; i++) {
		if (na[i] > 0) {
			double sum = 0.;
			l1 = l2;
			l2 += na[i];
			for (int l = l1; l < l2; l++) {
				sum += a[l] * x[ja[l]];
			}
			y[i] += sum;
		}
	}
}

// c--------------------------------------------------------------------------
// c  mode = 2  --  set x = x + a(transpose) * y
// c---------------------------------------------------------------------------
void mode2(int m, float *x, float *y, float *a, int *na, int *ja) {
	int l1 = 0, l2 = 0;
	for (int i = 0; i < m; i++) {
		if (na[i] > 0) {
			float yi = y[i];
			l1 = l2;
			l2 += na[i];
			for (int l = l1; l < l2; l++) {
				x[ja[l]] += (a[l] * yi);
			}
		}
	}
}
