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

void appendTime(Time *time_head, Time *new_time){
    Time *current = time_head;
    while(current != NULL){
        current = current->next;
    }
    current = new_time;
}

Event *createNewEvent(Earthquake eqk, char station_name_list[maxobs][MAXSTRLEN + 1], Time *obstime, int observationCnt){
    Event *new_event = (Event *)malloc(sizeof(Event));
	new_event->earthquake = eqk;
	memcpy(new_event->station_name_list, station_name_list, observationCnt * MAXSTRLEN+1);
    memcpy(new_event->observedTime, obstime, sizeof(Time) * observationCnt);
	new_event->next=NULL;
	return new_event;
}

Event *createEventList(char *filename){
    Time *obstimeList;
    FILE *fp_event = fopen(filename, "r");
    if(!fp_event){
		printf("Error on reading %s\n", filename);
		assert(0);
    }

    char str_eqk[100];
    if (fgets(str_eqk, sizeof(str_eqk), fp_event) == 0) {
        if (feof(fp_event)) {
            printf("eof on fp_event(%s)\n", filename);
            assert(0);
        } else {
            printf("error on reading fp_event(%s)\n", filename);
            assert(0);
        }
    }
    if(strlen(str_eqk) > sizeof(str_eqk)) {
        printf("str_inp is too long: %s\n", str_eqk);
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
    for(nsta = 0; fgets(str_eqk, sizeof(str_eqk), fp_event); nsta++) {
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
            printf("input length is too large. len=%d str_inp=%s\n", len,
                    str_eqk);
            assert(0);
        }

        char str_tmp[100];
        sscanf(str_eqk, "%s %d %d %d %d %f %99[^\n]\n", station_name[nsta], &iyr,
                &jday, &ihr, &imn, &sec, str_tmp);

        Time *time = createNewTime(iyr, jday, ihr, imn, sec);
        appendTime(obstimeList, time);

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
    fclose(fp_event);
}
