#include "sphrayderv/find.h"

void find(float *x, float xt, int n, int *ians) {
	int i1 = 0, im = 0, itst = 0, ihalf = 0;
	if (xt <= x[0]) goto a104;
	if (xt >= x[n]) goto a105;
	i1 = 1;
	im = n;
a102:
	itst = im - i1;
	if (itst > 1)
		goto a100;
	*ians = i1;
	goto a103;
a100:
	ihalf = (im + i1) / 2;
	if (xt > x[ihalf]) goto a101;
	im = ihalf;
	goto a102;
a101: i1 = ihalf;
a104:	*ians = 1;
	goto a103;
a105: *ians = n;
a103: return;
}
