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
    char phase[maxobs];
    float rwts[maxobs];
    char station_name_list[maxobs][MAXSTRLEN + 1];
}Event;

struct EventNode_{
    Event event;
    struct EventNode_ *next;
};
typedef struct EventNode_ EventNode;

TimeNode *createTimeNode(int, int, int, int, float);
void insertTime(TimeNode *, Time);
void appendTimeNode(TimeNode **, TimeNode *);
int getTimeCount(TimeNode *);
EventNode *createEventNode(Earthquake, char station_name_list[maxobs][MAXSTRLEN + 1], TimeNode *, char *, float *);
void appendEventNode(EventNode **, EventNode *);
int getEventCount(EventNode *);
EventNode *createEventList(char *);
Event *EventList2Arr(EventNode *);
int *checkTravelTime(Event , travelTimeTable *, StationNode *);
float *getObsTime(Event);
#endif