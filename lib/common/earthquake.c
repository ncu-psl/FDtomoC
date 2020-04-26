#include "common/earthquake.h"
Time *createNewTime(int iyr, int jday, int ihr, int imn, float sec){
    Time *new_time = (Time *)malloc(sizeof(Time));
    new_time->iyr = iyr;
    new_time->jday = jday;
    new_time->ihr = ihr;
    new_time->imn = imn;
    new_time->sec = sec;
    new_time->next = NULL;
    return new_time;
}

void appendTime(Time **time_head, Time *new_time){
    if(*time_head == NULL){
        *time_head = new_time;
        return;
    }
    Time *current = *time_head;
    while(current->next != NULL){
        current = current->next;
    }
    current->next = new_time;
}

Event *createNewEvent(Earthquake eqk, char station_name_list[maxobs][MAXSTRLEN + 1], Time *obstime, int observationCnt){
    Event *new_event = (Event *)malloc(sizeof(Event));
    new_event->observedTime = (Time *)malloc(sizeof(Time) * observationCnt);
	new_event->earthquake = eqk;
    new_event->observationCnt = observationCnt;
	memcpy(new_event->station_name_list, station_name_list, observationCnt * MAXSTRLEN+1);
    memcpy(new_event->observedTime, obstime, sizeof(Time) * observationCnt);
	new_event->next=NULL;
	return new_event;
}

void appendEvent(Event **event_head, Event *new_event){
    if(*event_head == NULL){
        *event_head = new_event;
        return;
    }
    Event *current =  *event_head;
    while (current->next != NULL){
        current = current->next;
    }
    current->next = new_event;
}


Event *createEventList(char *filename){
    Event *event_list;
    Event *event = NULL;

    FILE *fp_event = fopen(filename, "r");
    if(!fp_event){
		printf("Error on reading %s\n", filename);
		assert(0);
    }

    int eventCnt;
    for(eventCnt = 1; !feof(fp_event); eventCnt++){
        char str_eqk[100];
        if (fgets(str_eqk, sizeof(str_eqk), fp_event) == 0) {
            if(str_eqk[0] == '\n')
                break;
        }
        if(strlen(str_eqk) > sizeof(str_eqk)) {
            printf("str_eqk is too long: %s\n", str_eqk);
            assert(0);
        }

        int iyr, jday, ihr, imn;
        float sec;
        float xlat, xlon, dep;
        int eventId;

        sscanf(str_eqk, "%d %d %d %d %f %f %f %f %s\n", &iyr, &jday, &ihr,
                    &imn, &sec, &xlat, &xlon, &dep, eventId);

        Time time = {iyr, jday, ihr, imn, sec};
        Point3D location = {xlat, xlon, dep};
        Earthquake earthquake = {location, time};
        
        float pwt[maxobs], rwts[maxobs];
        char station_name[maxobs][MAXSTRLEN + 1];
        char phs[maxobs];
        int isgood[maxobs];
        char usemark;
        int nsta;
        Time *obstimeList;
        for(nsta = 0; fgets(str_eqk, sizeof(str_eqk), fp_event); nsta++) {
            if(str_eqk[0] == '\n')
                break;
            if (nsta >= maxobs) {
                printf("Error:  too many observations! nsta=%d maxobs=%d\n",
                        nsta, maxobs);
                assert(0);
            }
            
            trim(str_eqk);
            int len = strlen(str_eqk);
            if (str_eqk[0] == '\n') {
                break;
            }
            if (len > MAXSTRLEN) {
                printf("input length is too large. len=%d str_eqk=%s\n", len,
                        str_eqk);
                assert(0);
            }

            char str_tmp[100];
            sscanf(str_eqk, "%s %d %d %d %d %f %99[^\n]\n", station_name[nsta], &iyr,
                    &jday, &ihr, &imn, &sec, str_tmp);

            Time *time = createNewTime(iyr, jday, ihr, imn, sec);
            appendTime(&obstimeList, time);

            isgood[nsta] = 1;
            if (str_tmp[0] == '*') {
                isgood[nsta] = 0;
                usemark = '*';
                str_tmp[0] = ' ';
                trim(str_tmp);
            }
            sscanf(str_tmp, "%c %f", &phs[nsta], &rwts[nsta]);
    // c---convert from "0 1" to "P S" if necessary
            if (phs[nsta] == '0') {
                phs[nsta] = 'P';
            } else if (phs[nsta] == '1') {
                phs[nsta] = 'S';
            }
            pwt[nsta] = 1.f / (rwts[nsta] * rwts[nsta]);
        }
        event = createNewEvent(earthquake, station_name, obstimeList, nsta);
        appendEvent(&event_list, event);
    }
    fclose(fp_event);
    return event_list;
}

int *checkTravelTime(Event *event, travelTimeTable *table_list, Station *station_head){
    int numOfStations = getStationCount(station_head);
    int *timeIndex = (int *)malloc(sizeof(int) * event->observationCnt);
    travelTimeTable *current;
    for(int i = 0; i < event->observationCnt; i++){
        current = table_list;
        for (int j = 0; j < numOfStations; j++){
            if(strcmp(event->station_name_list[i], current->name) == 0){
                timeIndex[i] = j;
                break;
            }
            current++;
        }
    }
    return timeIndex;
}
