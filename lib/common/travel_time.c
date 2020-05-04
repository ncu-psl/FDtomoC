#include "common/travel_time.h"
void outputTravelTime(travelTimeTable time_table){
	FILE *fp_tmp = fopen(time_table.name, "wb");
	if(!fp_tmp){
		printf("Error happens while output travel time !\n");
		assert(0);
	}
	int size = sizeOfMesh3D(time_table.mesh);
	fwrite(time_table.time, sizeof(float), size, fp_tmp);
	fclose(fp_tmp);
}

