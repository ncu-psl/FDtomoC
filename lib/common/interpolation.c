#include "common/interpolation.h"

float linear_interpolation(float x, float a, float fa, float b, float fb){
        return (a == b) ? fa : fa + (fb-fa) * (x-a) / (b-a);
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
