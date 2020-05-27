#ifndef EARTHQUAKE
#define EARTHQUAKE
#include "common/grid.h"
#include "common/travel_time.h"
#include "common/station.h"
#include "common/time_process.h"

typedef struct{
    int iyr, jday, ihr, imn;
    float sec;
}Time; 

struct TimeNode_{
    Time time;
    struct TimeNode_ *next;
}; 

typedef struct TimeNode_ TimeNode;

typedef struct {
    Point3D location;
    Time time;
}Earthquake;

typedef struct {
    Earthquake earthquake;
    TimeNode *observedTimeList;
    char evid[maxobs];
    char phase[maxobs];
    int isgood[maxobs];
    float rwts[maxobs];
    char station_name_list[maxobs][MAXSTRLEN + 1];
}Event;

struct EventNode_{
    Event event;
    struct EventNode_ *next;
};
typedef struct EventNode_ EventNode;

typedef struct {
    float stdmin;
    double tmp_min[maxobs];
}LocData;

typedef struct {
    int iread, ivs, nthres, kmin, ndiv, ndiv2;
    float vpvs, resthres, resthrep, stdmax;
    char leqsfil[MAXSTRLEN + 1], fsumfil[MAXSTRLEN + 1], outlfil[MAXSTRLEN + 1],
                fhedfil[MAXSTRLEN + 1], fdatfil[MAXSTRLEN + 1];
}LocEnv;

TimeNode *createTimeNode(int, int, int, int, float);
void insertTime(TimeNode *, Time);
void appendTimeNode(TimeNode **, TimeNode *);
void copyTimeList(TimeNode **, TimeNode *);
TimeNode *TimeList2Arr(Time *, int);
int getTimeCount(TimeNode *);
EventNode *createEventNode(Earthquake, char station_name_list[maxobs][MAXSTRLEN + 1], TimeNode *,
                             char *, float *, int *, char *);
void appendEventNode(EventNode **, EventNode *);
int getEventCount(EventNode *);
EventNode *createEventList(char *);
Event *EventList2Arr(EventNode *);
int *checkTravelTime(Event , travelTimeTable *, int);
float *getObsTime(Event);
float *getPwt(Event);
void setLocFiles(LocEnv *, char *);
void setLocVariables(LocEnv *, char *);
LocEnv setLocEnv(char *);
#endif