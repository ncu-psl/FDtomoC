#ifndef STATION
#define STATION
#include "common/grid.h"
struct Station_{
    char *name;
    Point3D location;
};

struct StationNode_{
    int key;
    struct Station_ data;
    struct StationNode_ *next;
};

typedef struct Station_ Station;
typedef struct StationNode_ StationNode;

StationNode *createStationNode(char *, Point3D);
void insertStation(StationNode *, Station);
void appendStationNode(StationNode **, StationNode *);
StationNode *createStationList(char *, int);
int getStationCount(StationNode *);
#endif