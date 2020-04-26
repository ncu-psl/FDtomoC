#ifndef TRAVEL_TIME
#define TRAVEL_TIME
#include "common/grid.h"

typedef struct{
    char name[100];
    Mesh3D mesh;
    float *time;
}travelTimeTable;

#endif