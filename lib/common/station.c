#include "common/station.h"
Station *createNewStation(char *name, Point3D location){
    Station *new_station = (Station *)malloc(sizeof(Station));
	new_station->name = (char *)malloc(sizeof(char) * (strlen(name) + 1));
	new_station->location=location;
    memcpy(new_station->name, name, strlen(name) + 1);
	new_station->next=NULL;
	return new_station;
}

void appendStation(Station **station_head, Station *station) {
	if(*station_head == NULL){
		*station_head = station;
		return;
	}
	Station *ptr=*station_head;
	while(ptr->next!=NULL) {
		ptr=ptr->next;
	}
	ptr->next = station;
}


Station *createStationList(char *file, int sph){
    FILE *fp_sta = fopen(file, "r");
	if(!fp_sta) {
		printf("Can not open file: %s\n", file);
		assert(0);
	}
    Station *station_head = NULL;

    char str_inp[100];
	while(fgets(str_inp, sizeof(str_inp), fp_sta) != NULL){
        int zs = 0;
		float stx = 0, sty = 0, fxs = 0, fys = 0, fzs = 0;
		char sta[100];
		if(sph) {
			sscanf(str_inp, "%f %f %d %s %f %f", &sty, &stx, &zs, sta, &fys, &fxs);
		} else {
			sscanf(str_inp, "%f %f %f %s", &fys, &fxs, &fzs, sta);
		}
		fzs = zs / -1000.;
        
        Point3D location = {fxs, fys, fzs};
        Station *station =  createNewStation(sta, location);
        if(station_head==NULL){
            station_head = station;
        }else{
            appendStation(&station_head, station);
        }
    }
    return station_head;
}

int getStationCount(Station *station_list){
	int index = 0;
	Station *current = station_list;
	while(station_list != NULL){
		index++;
		station_list = station_list->next;
	}
	return index;
}
