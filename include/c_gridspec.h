//	common /gridspec/ x0, y0, z0, h, dq, df, x00, y00,
//     +        nxc, nyc, nzc, nx, ny, nz,
//     +        igridx(nxcm1),igridy(nycm1),igridz(nzcm1)
#include "c_parameter.h"
int nxc, nyc, nzc, nx, ny, nz;
double h, x0, y[1], z0, dq, df, x00, y00;
int igridx[nxcm1], igridy[nycm1], igridz[nzcm1];
