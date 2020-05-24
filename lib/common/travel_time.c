#include "common/travel_time.h"
void outputTravelTime(travelTimeTable time_table, char *filename){
	FILE *fp_tmp = fopen(filename, "wb");
	if(!fp_tmp){
		printf("Error happens while output travel time !\n");
		assert(0);
	}
	int size = sizeOfMesh3D(time_table.mesh);
	fwrite(time_table.time, sizeof(float), size, fp_tmp);
	fclose(fp_tmp);
}

void appendTableNode(travelTimeTableNode **table_head, travelTimeTableNode *new_table_node){
	if(*table_head == NULL){
        *table_head = new_table_node;
        return;
    }
    travelTimeTableNode *current =  *table_head;
    while (current->next != NULL){
        current = current->next;
    }
    current->next = new_table_node;
}

travelTimeTable createTable(Mesh3D mesh, char *name, float *time){
	travelTimeTable table;
	int size = sizeOfMesh3D(mesh);
	table.time = malloc(sizeof(float) * size);
	copyMesh3D(&table.mesh, &mesh);
	memcpy(table.name, name, strlen(name) + 1);
	memcpy(table.time, time, sizeof(float) * size);
	return table;
}
