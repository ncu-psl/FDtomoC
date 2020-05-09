#include "common/earthquake.h"
TimeNode *createTimeNode(int iyr, int jday, int ihr, int imn, float sec){
    TimeNode *new_time_node = (TimeNode *)malloc(sizeof(TimeNode));
    new_time_node->time.iyr = iyr;
    new_time_node->time.jday = jday;
    new_time_node->time.ihr = ihr;
    new_time_node->time.imn = imn;
    new_time_node->time.sec = sec;
    new_time_node->next = NULL;
    return new_time_node;
}

void insertTime(TimeNode *time_head, Time time){
    TimeNode *current = time_head;
	while(current->next!=NULL) {
		current=current->next;
	}
	current->time = time;
	current->next = NULL;
}


void appendTimeNode(TimeNode **time_head, TimeNode *new_time_node){
    if(*time_head == NULL){
        *time_head = new_time_node;
        return;
    }
    TimeNode *current = *time_head;
    while(current->next != NULL){
        current = current->next;
    }
    current->next = new_time_node;
}

int getTimeCount(TimeNode *time_list){
    if (time_list == NULL)  
        return 0;
	int index = 1;
	TimeNode *current = time_list;
	while(time_list->next != NULL){
		index++;
		time_list = time_list->next;
	}
	return index;
}

EventNode *createEventNode(Earthquake eqk, char station_name_list[maxobs][MAXSTRLEN + 1], TimeNode *obstime, 
                            char *phase, float *rwts, int *isgood, char *evid){
    EventNode *new_event_node = (EventNode *)malloc(sizeof(EventNode));
    int observationCnt = getTimeCount(obstime);
    if (observationCnt != 0){
        memcpy(new_event_node->event.station_name_list, station_name_list, observationCnt * (MAXSTRLEN+1));
        memcpy(new_event_node->event.phase, phase, observationCnt);
        memcpy(new_event_node->event.rwts, rwts, observationCnt * sizeof(float));
        memcpy(new_event_node->event.isgood, isgood, observationCnt * sizeof(int));
        memcpy(new_event_node->event.evid, evid, strlen(evid) + 1);
        new_event_node->event.observedTimeList = obstime;
    }else {
        new_event_node->event.observedTimeList = NULL;
    }
	new_event_node->event.earthquake = eqk;
	new_event_node->next=NULL;
	return new_event_node;
}

void appendEventNode(EventNode **event_head, EventNode *new_event_node){
    if(*event_head == NULL){
        *event_head = new_event_node;
        return;
    }
    EventNode *current =  *event_head;
    while (current->next != NULL){
        current = current->next;
    }
    current->next = new_event_node;
}

int getEventCount(EventNode *event_list){
    if (event_list == NULL)  
        return 0;
	int index = 1;
	EventNode *current = event_list;
	while(event_list->next != NULL){
		index++;
		event_list = event_list->next;
	}
	return index;
}

EventNode *createEventList(char *filename){
    EventNode *event_head = NULL;
    FILE *fp_event = fopen(filename, "r");
    if(!fp_event){
		printf("Error on reading %s\n", filename);
		assert(0);
    }

    int eventCnt;
    for(eventCnt = 1; !feof(fp_event); eventCnt++){
        char str_eqk[MAXSTRLEN];
        if (fgets(str_eqk, sizeof(str_eqk), fp_event) == 0) {
            break;
        }else if(str_eqk[0] == '\n'){
            break;
        }

        if(strlen(str_eqk) > sizeof(str_eqk)) {
            printf("str_eqk is too long: %s\n", str_eqk);
            assert(0);
        }

        int iyr, jday, ihr, imn;
        float sec;
        float xlat, xlon, dep;
        char  evid[11];

        sscanf(str_eqk, "%d %d %d %d %f %f %f %f %s\n", &iyr, &jday, &ihr,
                    &imn, &sec, &xlat, &xlon, &dep, evid);

        Time time = {iyr, jday, ihr, imn, sec};
        Point3D location = {xlat, xlon, dep};
        Earthquake earthquake = {location, time};
        
        float pwt[maxobs], rwts[maxobs];
        char station_name[maxobs][MAXSTRLEN + 1];
        char phs[maxobs];
        int isgood[maxobs];
        char usemark;
        int nsta;
        TimeNode *obstimeList = NULL;
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

            char str_tmp[MAXSTRLEN];
            sscanf(str_eqk, "%s %d %d %d %d %f %[^\n]\n", station_name[nsta], &iyr,
                    &jday, &ihr, &imn, &sec, str_tmp);

            TimeNode *new_time_node = createTimeNode(iyr, jday, ihr, imn, sec);
            appendTimeNode(&obstimeList, new_time_node);

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
        EventNode *new_event_node = createEventNode(earthquake, station_name, obstimeList, 
                                                    phs, rwts, isgood, evid);
        appendEventNode(&event_head, new_event_node);
    }
    fclose(fp_event);
    return event_head;
}

Event *EventList2Arr(EventNode *event_list){
    int event_count = getEventCount(event_list);
    Event *event_array = malloc(sizeof(Event) * event_count);
    for(int i = 0; i < event_count; i++){
        event_array[i] = event_list->event;
        event_array[i].observedTimeList = event_list->event.observedTimeList;
        event_list = event_list->next;
    }
    return event_array;
}



int *checkTravelTime(Event event, travelTimeTable *table_array, int table_size){
    int numbOfObservation = getTimeCount(event.observedTimeList);
    int *timeIndex = (int *)calloc(numbOfObservation, sizeof(int));
    travelTimeTable *current = NULL;
    for(int i = 0; i < numbOfObservation; i++){
        for (int j = 0; j < table_size; j++){
            if(strcmp(event.station_name_list[i], table_array[j].name) == 0){
                timeIndex[i] = j;
                break;
            }
        }
    }
    return timeIndex;
}

float *getObsTime(Event event){
    int numbOfObservation = getTimeCount(event.observedTimeList);
    float *obs_travel_time = calloc(numbOfObservation, sizeof(float));
    double origin_time = htoe2(event.earthquake.time);
    TimeNode *current_time = event.observedTimeList;
    for(int i = 0; i < numbOfObservation; i++){
        double obs_time = htoe2(current_time->time);
        obs_travel_time[i] = obs_time - origin_time;
        current_time = current_time->next;
    }
    return obs_travel_time;
}

float *getPwt(Event event){
    int numbOfObservation = getTimeCount(event.observedTimeList);
    float *pwt = malloc(sizeof(float) * numbOfObservation);
    for(int i = 0; i < numbOfObservation; i++){
        pwt[i] = 1.f / (event.rwts[i] * event.rwts[i]);
    }
    return pwt;
}

void setLocFiles(LocEnv *loc_env, char *spec_file){
    char *files[5] = {"leqsfil\0", "fsumfil\0", "outlfil\0", "fhedfil\0", "fdatfil\0" };
    char *loc_files[5] = {loc_env->leqsfil, loc_env->fsumfil, loc_env->outlfil, 
                            loc_env->fhedfil, loc_env->fdatfil};

    FILE *fp_spc;
    fp_spc = fopen(spec_file, "r");
	if (!fp_spc) {
		printf("(Error in read_spec.c)read fp_spc file error.\n");
		assert(0);
	}

    int len, ierr;
    char pval[MAXSTRLEN + 1];
    for (int i = 0; i < 5; i++) {
		get_vars(fp_spc, files[i], pval, &len, &ierr);
		if (ierr == 1) {
			printf("Error trying to read filename %s", files[i]);
			assert(0);
		}
		sscanf(pval, "%s", loc_files[i]);
	}

}

void setLocVariables(LocEnv *loc_env, char *spec_file){
    FILE *fp_spc;
    fp_spc = fopen(spec_file, "r");
	if (!fp_spc) {
		printf("(Error in read_spec.c)read fp_spc file error.\n");
		assert(0);
	}
    
    int len, ierr;
    char pval[MAXSTRLEN + 1];

    get_vars(fp_spc, "iread ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%d", &loc_env->iread);
	}
	if (loc_env->iread != 0 && loc_env->iread != 1) {
		loc_env->iread = 0;
	}
	get_vars(fp_spc, "ivs ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%d", &loc_env->ivs);
	}
	if (loc_env->ivs != 0 && loc_env->ivs != 1) {
		loc_env->ivs = 0;
	}
	get_vars(fp_spc, "vpvs ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &loc_env->vpvs);

	get_vars(fp_spc, "nthres ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &loc_env->nthres);
	get_vars(fp_spc, "resthres ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &loc_env->resthres);
	get_vars(fp_spc, "resthrep ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &loc_env->resthrep);
	get_vars(fp_spc, "stdmax ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%f", &loc_env->stdmax);
	get_vars(fp_spc, "kmin ", pval, &len, &ierr);
	if (ierr == 0) {
		sscanf(pval, "%d", &loc_env->kmin);
		loc_env->kmin--;
	}

//-----grid search control
	get_vars(fp_spc, "ndiv ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &loc_env->ndiv);
	if (loc_env->ndiv <= 0)
		loc_env->ndiv = 1;
	get_vars(fp_spc, "ndiv2 ", pval, &len, &ierr);
	if (ierr == 0)
		sscanf(pval, "%d", &loc_env->ndiv2);
	if (loc_env->ndiv2 <= 0)
		loc_env->ndiv2 = 1;
}

LocEnv setLocEnv(char *spec_file){
    LocEnv loc_env;
    loc_env.iread = 0;
//----ivs = 1 to treat Vp and Vs separately (individual time files).
//        = 0 to compute Ts as Tp*vpvs
	loc_env.ivs = 1;
	loc_env.vpvs = 1.78;
//---thresholds to define an acceptable location
//---number of phases threshold
	loc_env.nthres = 8;
//---residual threshold (absolute time)
	loc_env.resthres = .5;
//---residual threshold (percentage of travel time)
	loc_env.resthrep = 5.0;
//---std threshold
	loc_env.stdmax = 15.0;
	loc_env.kmin = 2;
    loc_env.ndiv = 20;
	loc_env.ndiv2 = 20;
    setLocFiles(&loc_env, spec_file);
    setLocVariables(&loc_env, spec_file);
    return loc_env;
}
