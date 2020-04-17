#ifndef GRID_H_
#define GRID_H_
#include "common/read_spec.h"

typedef struct {
    float x, y, z;
}Point3D; 

typedef struct {
    Point3D point;
    float value;
}Cell;

typedef struct{
    Point3D origin;
    int numberOfx,numberOfy, numberOfz;
    float xSpace, ySpace, zSpace;
    int *igridx, *igridy, *igridz;
}Mesh; 

int getNumberOfXfine(Mesh);
int getNumberOfYfine(Mesh);
int getNumberOfZfine(Mesh);
int sizeofFine(Mesh);
float *getXFineMesh(Mesh);
float *getYFineMesh(Mesh);
float *getZFineMesh(Mesh);
Point3D getFinePoint(Point3D,Mesh);
float *getXMesh(Mesh);
float *getYMesh(Mesh);
float *getZMesh(Mesh);
Point3D getCoarsePoint(Point3D point, Mesh mesh);
Mesh createMesh(SPEC spec);
void changeCoordinate(Mesh, int);
Point3D searchFineBase(Point3D, Mesh);
#endif