#ifndef STATION
#define STATION
#include "common/grid.h"
struct Station_{
    char *name;
    Point3DDouble location;
};

struct StationNode_{
    int key;
    struct Station_ data;
    struct StationNode_ *next;
};

typedef struct Station_ Station;
typedef struct StationNode_ StationNode;

StationNode *createStationNode(char *, Point3DDouble);
void insertStation(StationNode *, Station);
void appendStationNode(StationNode **, StationNode *);
StationNode *createStationList(char *, int);
int getStationCount(StationNode *);
Station *StationList2Arr(StationNode *);
#endif