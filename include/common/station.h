#ifndef STATION
#define STATION
#include "common/grid.h"
struct station{
    char *name;
    Point3D location;
    struct station *next;
};
typedef struct station Station;

Station *createNewStation(char *, Point3D);
void appendStation(Station **, Station *);
Station *createStationList(char *, int);
int getStationCount(Station *);
#endif