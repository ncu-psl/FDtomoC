#ifndef TIME_PROCESS
#define TIME_PROCESS
#include "common/earthquake.h"
#include <math.h>

int isleap(int);
void etoh(double, int *, int *, int *, int *, double *);
void htoe(int, int, int, int, double, double *);
void dtoepoch(int, double *);
double htoe2(Time);
#endif // TIME_PROCESS
