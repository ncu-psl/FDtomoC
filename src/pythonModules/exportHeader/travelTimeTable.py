#!/usr/bin/env python
from cffi import FFI
def travel_time_header():
    header = """
    typedef struct{
        char name[100];
        Mesh3D mesh;
        float *time;
    }travelTimeTable;
    """

    func = """
    travelTimeTable *sphfd(velocityModel3D, StationNode *);
    travelTimeTable sphfd_exec(velocityModel3D, Station);
    void outputTravelTime(travelTimeTable, char *);
    travelTimeTable createTable(Mesh3D, char *, float *);
    """

    return header + func