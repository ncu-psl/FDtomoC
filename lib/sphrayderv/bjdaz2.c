#include "sphrayderv/bjdaz2.h"

void bjdaz2(float zlat0, float alon0, float zlat, float alon, float* delt, float* az1, float* az2, int lattype) {
	double alat0, alat, st0, ct0, ct1, s0c1, st1, s1c0, dlon, sdlon, cdlon, cdelt, a, b, sdelt, ddelt;
	double aze = 0, twopi = pii * 2, azs = pii;
	if (lattype == 1) {
		alat0 = pii * 0.5 - zlat0;
		alat = pii * 0.5 - zlat;
	}
	else {
		alat0 = zlat0;
		alat = zlat;
	}
	st0 = sin(alat0);
	ct0 = cos(alat0);
	ct1 = cos(alat);
	s0c1 = st0 * ct1;
	st1 = sin(alat);
	s1c0 = st1 * ct0;
	dlon = alon - alon0;
	sdlon = sin(dlon);
	cdlon = cos(dlon);
	cdelt = st0 * st1*cdlon + ct0 * ct1;
	b = s0c1 - s1c0 * cdlon;
	a = st1 * sdlon;
	sdelt = sqrt(b*b + a * a);
	ddelt = atan2(sdelt, cdelt);
	if (sdelt != 0) aze = atan2(a, b);
	a = -sdlon * st0;
	b = s1c0 - s0c1 * cdlon;
	azs = pii;
	if (sdelt != 0) azs = atan2(a, b);
	if (aze < 0) aze = aze + twopi;
	if (azs < 0) azs = azs + twopi;
	*delt = ddelt;
	*az1 = aze;
	*az2 = azs;
	return;
}