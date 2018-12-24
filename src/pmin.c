#include "../include/sphrayderv/pmin.h"

double pmin(double d1, double d2) {
	// --- return minimum of positive d1 and d2.If neither are positive return -1e20
	if (d1 >= 0) {
		if (d2 >= 0) {
			if (d1 < d2) {
				return d1;
			}
			else {
				return d2;
			}
		}
		else {
			return d1;
		}
	}
	else {
		if (d2 >= 0) {
			return d2;
		}
		else {
			return -1e20;
		}
	}
}