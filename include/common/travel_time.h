#ifndef TRAVEL_TIME
#define TRAVEL_TIME
#include "common/grid.h"
#include "common/parameter.h"

typedef struct{
    char name[100];
    Mesh3D mesh;
    float *time;
}travelTimeTable;

typedef struct travelTimeTableNode_{
    travelTimeTable table;
    struct travelTimeTableNode_ *next;
};

typedef struct travelTimeTableNode_ travelTimeTableNode;
void appendTableNode(travelTimeTableNode **, travelTimeTableNode *);
void outputTravelTime(travelTimeTable, char *);
travelTimeTable createTable(Mesh3D, char *, float *);
#endif