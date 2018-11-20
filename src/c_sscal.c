#include "../include/c_sscal.h"

// replace single precision sx by single precision sa*sx.
// for i = 0 to n-1, replace sx(1+i*incx) with  sa * sx(1+i*incx)
void sscal(int n, float sa, float *sx, int incx) {
	if (n <= 0)
		return;
	if (incx == 1) {
// c        code for increments equal to 1.
// c        clean-up loop so remaining vector length is a multiple of 5.
		int m = n % 5;
		if (m != 0) {
			for (int i = 0; i < m; i++) {
				sx[i] *= sa;
			}
			if (n < 5) {
				return;
			}
		}
		int mp1 = m;
		for (int i = mp1; i < n; i += 5) {
			for (int j = 0; j < 5; j++) {
				sx[i + j] *= sa;
			}
		}
	} else {
// ode for increments not equal to 1.
		int ns = n * incx;
		for (int i = 0; i < ns; i += incx) {
			sx[i] *= sa;
		}
	}
}
