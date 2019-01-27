#include "runlsqr/normlz.h"
#include "runlsqr/snrm2.h"
#include "runlsqr/sscal.h"

//      normlz  is required by subroutine lsqr.  it computes the
//      euclidean norm of  x  and returns the value in  beta.
//      if  x  is nonzero, it is scaled so that norm(x) = 1.
//
//      functions and subroutines
//
//      blas       snrm2,sscal
void normlz(int n, float *x, float *beta) {
	*beta = snrm2(n, x, 1);
	if (*beta > 0) {
		sscal(n, 1 / *beta, x, 1);
	}
}
