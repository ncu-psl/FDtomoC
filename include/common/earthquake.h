#ifndef EARTHQUAKE
#define EARTHQUAKE
#include "common/grid.h"
#include "common/travel_time.h"
#include "common/station.h"
struct time{
    int iyr, jday, ihr, imn;
    float sec;
    struct time *next;
}; 
typedef struct time Time;

typedef struct {
    Point3D location;
    Time time;
}Earthquake;

struct event{
    Earthquake earthquake;
    Time *observedTime;
    char station_name_list[maxobs][MAXSTRLEN + 1];
    int eventId, observationCnt;
    struct event *next; 
};
typedef struct event Event;

Time *createNewTime(int, int, int, int, float);
void appendTime(Time **, Time *);
Event *createNewEvent(Earthquake, char station_name_list[maxobs][MAXSTRLEN + 1], Time *, int);
void appendEvent(Event **, Event *);
Event *createEventList(char *);
int *checkTravelTime(Event *, travelTimeTable *, Station *);
#endif