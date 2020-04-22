#ifndef GRID_H_
#define GRID_H_
#include "vec/vec.h"
#include "common/read_spec.h"

typedef struct {
    float x, y, z;
}Point3D; 

typedef struct {
    Point3D point;
    float value;
}Cell;

typedef struct{
    int numberOfNode;
    int space;
    vec_int_t igrid;
}Mesh1D; 

typedef struct{
    Point3D numberOfNode;
    Point3D space;
    int *igridx;
    int *igridy;
    int *igridz;
}Mesh3D;

typedef struct{
    Mesh1D mesh1d;
    int unit;
    int origin;
}Coordinate1D;

typedef struct{
    Mesh3D mesh;
    Point3D space;
    Point3D origin;
}Coordinate3D;

int sizeOfMesh3D(Mesh3D);
float *getXAxis(Coordinate3D);
float *getYAxis(Coordinate3D);
float *getZAxis(Coordinate3D);
Point3D getPoint3D(Point3D, Coordinate3D);
Mesh1D createMesh1D(int, int, int*);
Coordinate3D change2Sphere(Coordinate3D, int);
Point3D searchFineBase(Point3D, Coordinate3D);
#endif