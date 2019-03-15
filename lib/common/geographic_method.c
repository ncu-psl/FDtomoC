#include "common/geographic_method.h"

// c------convert geographic latitude to geocentric latitude-------------
// c        hlat (input) = geographic latitude in radians (north positive)
// c        glat (output)= geocentric latitude in radians (north positive)
// c---------------------------------------------------------------------
double glat(double hlat) {
	double halfpi = 1.570796f, polfac = 0.010632f, elfac = 0.993277f;
	if (halfpi - fabs(hlat) >= (double) 0.05) {
		return atan(elfac * sin(hlat) / cos(hlat));
	}
	else {
		//c------special formula near pole
		return hlat / elfac - fabs(polfac) * (hlat >= 0 ? 1 : -1);
	}
}

// c------convert geocentric latitude back to geographic latitude------------ -
// c       glatinv(output) = geographic latitude in radians(north positive)
// c       hlat(input) = geocentric latitude in radians(north positive)
// c-------------------------------------------------------------------- -
double glatinv(double hlat) {
	double halfpi = 1.570796;
	if (halfpi - fabs(hlat) >= 0.05) {
		return atan(sin(hlat) / cos(hlat) / 0.993277);
	}
	else {
		if (hlat > 0) {
			return (hlat + 0.010632)*0.993277;
		}
		else {
			return (hlat - 0.010632)*0.993277;
		}
	}
}

// c--------------------------------------------------------------------
// c
// c       glath.f
// c
// c       This is a version of glat that takes into account
// c       changes in elevation.
// c
// c       This version uses the WGS84 ellipsoid
// c
// c       Input
// c               xlat    Geodetic/Geogrphic latitude in radians (north positive)
// c               h       Ellipsoidal depth in km.  Note this is the usual GPS height
// c                       rather than the often used Geoidal elevation (MSL).
// c
// c       Output
// c               glath   Geocentric latitude in radians (returned)
// c               r       Distance from center of the Earth in km
// c
// c       All operations are double precision
// c
// c       Author:  S. Roecker, RPI.   July, 2008
// c

double glath(double xlat, double h, double *r) {
	//c---ap is semi major axis, bp is semiminor axis, f is inverse flattening, esq is the square of
	//c       the ellipticity, which we compute from fi*(2.d0-fi) where fi is 1/f.

	double ap = 6378137.0; //, bp = 6356752.314245, f = 298.257223563;
	double esq = 6.69437978616733379 * 0.001;	//-03;

	double sinxl = sin(xlat);

	//c---convert depth in km to elevation in meters
	double hm = -h * 1000.0;

	//c       anu is the ellipsoidal radius of curvature at the current geographic latitude
	double anu = ap / sqrt(1.0 - esq * sinxl * sinxl);

	double x = (anu + hm) * cos(xlat);
	double z = ((1.0 - esq) * anu + hm) * sinxl;

	*r = sqrt(x * x + z * z) / 1000.0;
	return atan2(z, x);

}

//c--------------------------------------------------------------------
//c
//c       glathinv.f
//c
//c
//c       This is a version of glatinv that takes into account
//c       changes in elevation.  This does the reverse of glath.
//c
//c       This version uses the WGS84 ellipsoid
//c
//c       Input
//c               xcent   Geocentric Latitude in radians
//c               r       Distance from center of Earth in km
//c
//c       Output
//c               glathinv Geodetic/Geogrphic latitude in radians (returned)
//c               h       Ellipsoidal depth in km.  Note this is the usual GPS height
//c                       rather than the often used Geoidal elevation (MSL).
//c
//c       All operations are double precision
//c
//c       Author:  S. Roecker, RPI.   July, 2008
//c

double glathinv(double xcent, double r, double *h) {

	//implicit real*8 (a-h,o-z)

	//c---ap is semi major axis, bp is semiminor axis, f is inverse flattening, esq is the square of
	//c       the ellipticity, which we compute from fi*(2.d0-fi) where fi is 1/f.
	double ap = 6378137.0;    //, bp = 6356752.314245, f = 298.257223563;
	double esq = 6.69437978616733379 * 0.001;
	//c---After about 5 iterations, the precision should be on the order of cm. We should
	//c       probably set a tolerance for this in any event
	int maxitt = 5;

	//c---convert to meters
	double rm = r * 1000.0;

	double xlatr = xcent;
	double x = rm * cos(xlatr);
	double z = rm * sin(xlatr);

	double anu;

	for (int nitt = 0; nitt < maxitt; nitt++) {
		double sinxl = sin(xlatr);
		anu = ap / sqrt(1.0 - esq * sinxl * sinxl);
		double hb = x / cos(xlatr) - anu;
		double zp = z / (1.0 - esq * (anu / (anu + hb)));
		xlatr = atan2(zp, x);
		//c         write(*,*) xlatr/degrad, hb
	}
	*h = x / cos(xlatr) - anu;
	*h = -*h / 1000.0;
	return xlatr;
}
