#!/usr/bin/env python
from cffi import FFI
def event_header():
    header = """
    #define SIZEOFA 200000000
    #define NMAX 2500000
    #define maxobs 800
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

    typedef struct{
        float elements[SIZEOFA];
        int column_elements[SIZEOFA];
        int elements_row[NMAX];
        int *jndx;
        int number_columns;
        int number_rows;
        int total_elements;
    }sparse_matrix;

    typedef struct{
        sparse_matrix *mat;
        float *b;
    }SPHRAYDERV_DATA;

    """

    func = """
    TimeNode *createTimeNode(int, int, int, int, float);
    void insertTime(TimeNode *, Time);
    void appendTimeNode(TimeNode **, TimeNode *);
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
    Event singleLoc(Coordinate3D, travelTimeTable *, Event, int, LocEnv);
    EventNode *sphfdloc(Coordinate3D, travelTimeTable *, int, Event *, int, LocEnv);
    SPHRAYDERV_DATA *sphrayderv(velocityModel3D, travelTimeTable *, Event *, int, Station *, int, SphraydervEnv, CommonEnv);
    RUNLSQR_DATA *runlsqr(SPHRAYDERV_DATA *, RunlsqrEnv, CommonEnv);
    """

    return header + func
