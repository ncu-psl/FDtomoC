#include "common/station.h"
StationNode *createStationNode(char *name, Point3DDouble location){
    StationNode *new_station_node = (StationNode *)malloc(sizeof(StationNode));
	new_station_node->data.name = (char *)malloc(sizeof(char) * (strlen(name) + 1));
	new_station_node->data.location = location;
    memcpy(new_station_node->data.name, name, strlen(name) + 1);
	new_station_node->next=NULL;
	return new_station_node;
}

void insertStation(StationNode *station_head, Station station) {
	StationNode *current = station_head;
	while(current->next!=NULL) {
		current=current->next;
	}
	current->data = station;
	current->next = NULL;
}

void appendStationNode(StationNode **station_head, StationNode *new_station_node) {
	if(*station_head == NULL){
		*station_head = new_station_node;
		return;
	}
	StationNode *current = *station_head;
	while(current->next != NULL) {
		current = current->next;
	}
	current->next = new_station_node;
}


StationNode *createStationList(char *file, int sph){
    FILE *fp_sta = fopen(file, "r");
	if(!fp_sta) {
		printf("Can not open file: %s\n", file);
		assert(0);
	}
	StationNode *station_head = NULL;
    char str_inp[100];
	while(fgets(str_inp, sizeof(str_inp), fp_sta) != NULL){
        int zs = 0;
		float stx = 0, sty = 0;
		double fxs = 0, fys = 0, fzs = 0;
		char sta[100];
		if(sph) {
			sscanf(str_inp, "%f %f %d %s %lf %lf", &sty, &stx, &zs, sta, &fys, &fxs);
		} else {
			sscanf(str_inp, "%lf %lf %lf %s", &fys, &fxs, &fzs, sta);
		}
		fzs = zs / -1000.;
        
		Point3DDouble location = {fxs, fys, fzs};
        StationNode *new_station_node = createStationNode(sta, location);
    	appendStationNode(&station_head, new_station_node);
    }
    return station_head;
}

int getStationCount(StationNode *station_list){
	if (station_list == NULL)
		return 0;
	int index = 1;
	StationNode *current = station_list;
	while(station_list->next != NULL){
		index++;
		station_list = station_list->next;
	}
	return index;
}

Station *StationList2Arr(StationNode *station_list){
	int station_count = getStationCount(station_list);
    Station *station_array = malloc(sizeof(Station) * station_count);
	StationNode *current_station = station_list;
    for(int i = 0; i < station_count; i++){
        station_array[i] = current_station->data;
        current_station = current_station->next;
    }
    return station_array;
}
