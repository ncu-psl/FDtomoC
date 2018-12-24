#include "../include/sphrayderv/elpcr.h"
#include "../include/sphrayderv/find.h"
void elpcr(double elat, float delta, float azdg, float depth, float *elpc, int iphase, int ic) {
	// c    A set of subroutines useful for spherical - earth calculations
	// c    ellipticity subs modified from ken creager's subs
	// c
	// c
	// c   Compute ellipticity correction from dziewonski and gilbert
	// c   Geophys.J.R.Astr.Soc. (1976) 44, 7 - 17.
	// c   elat = epicentral latitude
	// c   delta = distance in degrees
	// c   azdg = azimuth from epicenter to receiver in degrees
	// c   depth = earthquake depth in km.
	// c   iphase = 0 for p or 1 for s
	// c   elpc = ellipticity correction in seconds
	// c   if input latitude is in geographic(? geocentric ? ) degrees   ic = 0
	// c   if input latitude is in geocentric radians   ic = 1	

	if (delta < 0 || delta>180) {
		goto a999;
	}

	// interpolate over distance and depth to get tau0, tau1, and tau2.

	int idl = 0;
	double y = 0;
	if (iphase == 1) {
		if (delta > 105) {
			goto a999;
		}
		find(ecsd, delta, 22, &idl);
		y = (delta - ecsd[idl]) / (ecsd[idl + 1] - ecsd[idl]);
	}
	else {
		if (delta > 180) {
			goto a999;
		}
		find(ecpd, delta, 37, &idl);
		y = (delta - ecpd[idl]) / (ecpd[idl + 1] - ecpd[idl]);
	}
	int idt = 0;
	double x = 0;
	if (depth <= 300) {
		idt = 0;
		x = depth / 300;
		if (x < 0) {
			x = 0;
		}
	}
	else {
		idt = 1;
		x = (depth - 300) / 350;
		if (x > 1) {
			x = 1;
		}
	}
	double a = 1 - x - y + x * y;
	double b = x - x * y;
	double c = y - x * y;
	double d = x * y;

	double tau[3];
	if (iphase == 1) {
		tau[0] = ecs[0][idl][idt] * a + ecs[0][idl][idt + 1] * b + ecs[0][idl + 1][idt] * c + ecs[0][idl + 1][idt + 1] * d;
		tau[1] = ecs[1][idl][idt] * a + ecs[1][idl][idt + 1] * b + ecs[1][idl + 1][idt] * c + ecs[1][idl + 1][idt + 1] * d;
		tau[2] = ecs[2][idl][idt] * a + ecs[2][idl][idt + 1] * b + ecs[2][idl + 1][idt] * c + ecs[2][idl + 1][idt + 1] * d;
	}
	else {
		tau[0] = ecp[0][idl][idt] * a + ecp[0][idl][idt + 1] * b + ecp[0][idl + 1][idt] * c + ecp[0][idl + 1][idt + 1] * d;
		tau[1] = ecs[1][idl][idt] * a + ecs[1][idl][idt + 1] * b + ecs[1][idl + 1][idt] * c + ecs[1][idl + 1][idt + 1] * d;
		tau[2] = ecs[2][idl][idt] * a + ecs[2][idl][idt + 1] * b + ecs[2][idl + 1][idt] * c + ecs[2][idl + 1][idt + 1] * d;
	}

	double theta = elat;
	if (ic == 0) {
		theta *= PI / 180;
	}
	theta = PI / 2 - theta;
	double az = azdg * PI / 180;
	double cc[3];
	cc[0] = (1 + 3 * cos(2 * theta)) / 4;
	cc[1] = sin(2 * theta)*cos(az)*sqrt(3) / 2;
	cc[2] = sin(theta)*sin(theta)*cos(2 * az)*sqrt(3) / 2;
	*elpc = 0;
	for (int i = 0; i < 3; i++) {
		*elpc += cc[i] * tau[i];
	}
a999:
	*elpc = 0;
	return;
}
