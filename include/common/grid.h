#ifndef GRID_H_
#define GRID_H_
#include "common/geographic_method.h"
#include "common/vec/vec.h"
#include "common/read_spec.h"

typedef struct {
    float x, y, z;
}Point3D; 

typedef struct {
    double x, y, z;
}Point3DDouble;

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
    double xspace;
    double yspace;
    double zspace;
    vec_int_t igridx;
    vec_int_t igridy;
    vec_int_t igridz;
}Mesh3D;

typedef struct{
    Mesh1D mesh1d;
    int unit;
    int origin;
}Coordinate1D;

typedef struct{
    Mesh3D mesh;
    Point3DDouble space;
    Point3DDouble origin;
}Coordinate3D;

int sizeOfMesh3D(Mesh3D);
float *getAxis(Coordinate1D);
float *getXAxis(Coordinate3D);
float *getYAxis(Coordinate3D);
float *getZAxis(Coordinate3D);
Point3D getPoint3D(Point3D, Coordinate3D);
Mesh1D createMesh1D(int, int, vec_int_t);
Mesh3D readFineMesh3D(SPEC);
Mesh3D readCoarseMesh3D(SPEC);
Coordinate1D createCoordinate(Mesh1D, int, int);
Coordinate3D readFineCoordinate(SPEC);
Coordinate3D readCoarseCoordinate(SPEC);
Coordinate3D change2Sphere(Coordinate3D, int);
Point3D searchFineBase(Point3D, Coordinate3D);
#endif