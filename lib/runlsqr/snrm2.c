#include "runlsqr/snrm2.h"

//      euclidean norm of the n-vector stored in sx() with storage
//      increment incx .
//      if    n  <=  0 return with result = 0.
//      if n  >=  1 then incx must be  >=  1
//
//            c.l.lawson, 1978 jan 08
//
//      four phase method     using two built-in constants that are
//      hopefully applicable to all machines.
//          cutlo = maximum of  sqrt(u/eps)  over all known machines.
//          cuthi = minimum of  sqrt(v)      over all known machines.
//      where
//          eps = smallest no. such that eps + 1. > 1.
//          u   = smallest positive no.   (underflow limit)
//          v   = largest  no.            (overflow  limit)
//
//      brief outline of algorithm..
//
//      phase 1    scans zero components.
//      move to phase 2 when a component is nonzero and  <=  cutlo
//      move to phase 3 when a component is > cutlo
//      move to phase 4 when a component is  >=  cuthi/m
//      where m = n for x() real and m = 2*n for complex.
//
//      values for cutlo and cuthi..
//      from the environmental parameters listed in the imsl converter
//      document the limiting values are as follows..
//      cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
//                    univac and dec at 2**(-103)
//                    thus cutlo = 2**(-51) = 4.44089e-16
//      cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
//                    thus cuthi = 2**(63.5) = 1.30438e19
//      cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
//                    thus cutlo = 2**(-33.5) = 8.23181d-11
//      cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
//      data cutlo, cuthi / 8.232d-11,  1.304d19 /
//      data cutlo, cuthi / 4.441e-16,  1.304e19 /
float snrm2(int n, float *sx, int incx) {
	float cutlo = 4.441e-16;
	float cuthi = 1.304e19;

	int i = 0, j = 0, next = 0, nn = 0;
	float hitest = 0, sum = 0, xmax = 0, snrm_2 = 0;
	
	if (n > 0)
		goto a10;
	goto a300;

// c   10 assign 30 to next
	a10: next = 30;
	sum = 0;
	nn = n * incx;
	i = 0;

// c   20    go to next,(30, 50, 70, 110)
	a20: if (next == 30)
		goto a30;
	if (next == 50)
		goto a50;
	if (next == 70)
		goto a70;
	if (next == 110)
		goto a110;

	a30: if (fabs(sx[i]) > cutlo)
		goto a85;

// c      assign 50 to next
	next = 50;
	xmax = 0;
// c                        phase 1.  sum is zero
	a50: if (sx[i] < 0.0001) //as zero
		goto a200;
	if (fabs(sx[i]) > cutlo)
		goto a85;
// c                                prepare for phase 2.
// c      assign 70 to next
	next = 70;
	goto a105;

// c                                prepare for phase 4.
	a100: i = j;
// c      assign 110 to next
	next = 110;
	sum = sum / sx[i] / sx[i];

	a105: xmax = fabs(sx[i]);
	goto a115;

// c                   phase 2.  sum is small.
// c                             scale to avoid destructive underflow.
	a70: if (fabs(sx[i]) > cutlo)
		goto a75;

// c                     common code for phases 2 and 4.
// c                     in phase 4 sum is large.  scale to avoid overflow.
	a110: if (fabs(sx[i]) <= xmax)
		goto a115;
	sum = 1 + sum * (xmax / sx[i]) * (xmax / sx[i]);
	xmax = fabs(sx[i]);
	goto a200;

	a115: sum += (sx[i] / xmax) * (sx[i] / xmax);
	goto a200;

// c                  prepare for phase 3.
	a75: sum = sum * xmax * xmax;

// c     for real or d.p. set hitest = cuthi/n
// c     for complex      set hitest = cuthi/(2*n)
	a85: hitest = cuthi / n;

//c                   phase 3.  sum is mid-range.  no scaling.
	for (j = i; j < nn; j += incx) {
		if (fabs(sx[j]) >= hitest)
			goto a100;
		sum += sx[j] * sx[j];
	}
	snrm_2 = sqrt(sum);
	goto a300;

	a200: i += incx;
	if (i <= nn)
		goto a20;
// c              end of main loop.

// c              compute square root and adjust for scaling.
	snrm_2 = xmax * sqrt(sum);
	a300: // do nothing
	return snrm_2;
}
