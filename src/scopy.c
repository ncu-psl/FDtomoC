#include "../include/runlsqr/scopy.h"

//      copy single precision sx to single precision sy.
//      for i = 0 to n-1, copy  sx(lx+i*incx) to sy(ly+i*incy),
//      where lx = 1 if incx  >=  0, else lx = (-incx)*n, and ly is
//      defined in a similar way using incy.
void scopy(int n, float *sx, int incx, float *sy, int incy) {
	int m = 0;
	if (n <= 0)
		return;
	if (incx == incy) {
		if (incx < 1) {
			goto a5;
		} else if (incx == 1) {
			goto a20;
		} else {
			goto a60;
		}
	}

// c        code for unequal or nonpositive increments.

	a5: {
		int ix = 0, iy = 0;
		if (incx < 0)
			ix = (-n + 1) * incx;
		if (incy < 0)
			iy = (-n + 1) * incy;
		for (int i = 0; i < n; i++) {
			sy[iy] = sx[ix];
			ix += incx;
			iy += incy;
		}
		return;
	}

// c        code for both increments equal to 1
// c        clean-up loop so remaining vector length is a multiple of 7.

	a20: {
		m = n % 7;
		if (m == 0)
			goto a40;
		for (int i = 0; i < m; i++) {
			sy[i] = sx[i];
		}
		if (n < 7)
			return;

		a40: {
			int mp1 = m;
			for (int i = mp1; i < n; i += 7) {
				for (int j = 0; j < 7; j++) {
					sy[i + j] = sx[i + j];
				}
			}
			return;
		}
	}

// c        code for equal, positive, nonunit increments.

	a60: {
		int ns = n * incx;
		for (int i = 0; i < ns; i += incx) {
			sy[i] = sx[i];
		}
		return;
	}
}
