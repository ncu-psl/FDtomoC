#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include 	<assert.h>
#include    <string.h>
#include    <fcntl.h>
#define nhbyte 58 * 4

struct vhead {
	char header[120];
	double fxs, fys, fzs;
	double clat, clon, cz;
	double x0, y0, z0, dx, dy, dz;
	float az;
	int nx, ny, nz;
};
int main(int argc, char *argv[]){
	char fmod[100], cmod[100];
    struct vhead cheadin, fheadin;
	float *cvsave, *fvsave;
	int nxyz;

	strcpy(fmod, argv[1]);
	strcpy(cmod, argv[2]);

	sscanf(argv[3], "%d", &nxyz);
	int nxyz2 = nxyz*2;
	cvsave = (float *) malloc(sizeof(float) * nxyz2);
	fvsave = (float *) malloc(sizeof(float) * nxyz2);
   
	
    FILE *fp_cmod=fopen(cmod, "rb");
    if (!fp_cmod) {
		printf("Cannot open cmod: %s\n", cmod);
		assert(0);
	}
	fread(&cheadin, sizeof(char), nhbyte, fp_cmod);
    fread(cvsave, sizeof(cvsave[0]), nxyz2, fp_cmod);
	  
    FILE *fp_fmod=fopen(fmod, "rb");
    if (!fp_fmod) {
		printf("Cannot open cmod: %s\n", fmod);
		assert(0);
	}
	fread(&fheadin, sizeof(char), nhbyte, fp_fmod);
    fread(fvsave, sizeof(fvsave[0]), nxyz2, fp_fmod);
    
	/*
	FILE *fp_fnw = fopen("transaltion.txt", "w");
	fprintf(fp_fnw, "header : %s\n", cheadin.header);
	fprintf(fp_fnw, "fxs = %lf, fys = %lf, fzs = %lf\n", cheadin.fxs, cheadin.fys, cheadin.fzs);
	fprintf(fp_fnw, "clat = %lf, clon = %lf, cz = %lf\n", cheadin.clat, cheadin.clon, cheadin.cz);
	fprintf(fp_fnw, "x0 = %lf, y0 = %lf, z0 = %lf\n", cheadin.x0, cheadin.y0, cheadin.z0);
	fprintf(fp_fnw, "dx = %lf, dy = %lf, dz = %lf\n", cheadin.dx, cheadin.dy, cheadin.dz);
	*/
	FILE *fp_comp = fopen("compare.txt", "w");
	int flag = 0;
	for(int i = 0; i < nxyz2; i++){
		if (fvsave[i] !=  cvsave[i]){
			fprintf(fp_comp, "fmod[%d] = %f  ||||  cmod[%d] = %f\n", i, fvsave[i], i, cvsave[i]);
			flag = 1;
		}
	}
	if(flag)
		printf("fmod is different from cmod !!\n");
	fclose(fp_comp);

}
