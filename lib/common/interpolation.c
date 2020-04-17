#include "common/velocity_model.h"
#include "common/interpolation.h"
float linear_interpolation(float x, float a, float fa, float b, float fb){
        return (a == b) ? fa : fa + (fb-fa) * (x-a) / (b-a);
}

float trilinear_interpolation(Point3D point, Cell cells[2][2][2]){
	float x0 = cells[0][0][0].point.x;
	float x1 = cells[1][0][0].point.x;
	float y0 = cells[0][0][0].point.y;
	float y1 = cells[0][1][0].point.y;
	float z0 = cells[0][0][0].point.z;
	float z1 = cells[0][0][1].point.z;

	float xd = (point.x - x0) / x1 - x0;
	float yd = (point.y - y0) / y1 - y0;
	float zd = (point.z - z0) / z1 - z0;

	float c00 = cells[0][0][0].value * (1 - xd) + cells[1][0][0].value * xd;
	float c01 = cells[0][0][1].value * (1 - xd) + cells[1][0][1].value * xd;
	float c10 = cells[0][1][0].value * (1 - xd) + cells[1][1][0].value * xd;
	float c11 = cells[0][1][1].value * (1 - xd) + cells[1][1][1].value * xd;

	float c0 = c00 * (1 - yd) + c10 * yd;
	float c1 = c01 * (1 - yd) + c11 * yd;

	float c = c0 * (1 - yd) + c1 * zd;
	return c;
}

float *linear_interpolation_array(float *c, float *x, float *y, int n, int k, char *mode){
    float *fc = (float *)malloc(sizeof(float) * k);
    for(int i = 0;i < k; i++){
        int ik;
		for (ik = 1; ik < n; ik++) {
			if (x[ik] > c[i])
				break;
		}
		ik--;
		float v = 0;
        if (mode[i] == 'I') {
			fc[i] = linear_interpolation(c[i], x[ik], y[ik], x[ik + 1], y[ik + 1]);
		} else {
			fc[i] = y[ik];
		}
    }
    return fc;
}

float trilinear_interpolation_base(Point3D point, Point3D base, velocity3D model){
	base = getFinePoint(base, model.mesh);
	/*
    Cell c000 = {base,  model.};
    Cell c001 = {gx[i], gy[j], gz[k + 1]};
    Cell c010 = {gx[i], gy[j + 1], gz[k]};
    Cell c011 = {gx[i], gy[j + 1], gz[k + 1]};
    Cell c100 = {gx[i + 1], gy[j], gz[k]};
    Cell c101 = {gx[i + 1], gy[j], gz[k + 1]};
    Cell c110 = {gx[i + 1], gy[j + 1], gz[k]};
    Cell c111 = {gx[i + 1], gy[j + 1], gz[k + 1]};
    Cell cells[2][2][2] = {c000 , c001, c010, c011, c100, c101, c110, c111};
    float tmp = trilinear_interpolation(point, cells);
	*/
}
