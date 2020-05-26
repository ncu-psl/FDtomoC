#!/usr/bin/env python
from cffi import FFI
def grid_header():
    header = """
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
        int *igrid;
    }Mesh1D; 

    typedef struct{
        Point3D numberOfNode;
        int *gridx;
        int *gridy;
        int *gridz;
    }Mesh3D;

    typedef struct{
        Mesh1D mesh;
        int space;
        int origin;
    }Coordinate1D;

    typedef struct{
        Mesh3D mesh;
        Point3DDouble space;
        Point3DDouble origin;
    }Coordinate3D;
    """

    func = """
    int sizeOfMesh3D(Mesh3D);
    float *getAxis(Coordinate1D);
    float *getXAxis(Coordinate3D);
    float *getYAxis(Coordinate3D);
    float *getZAxis(Coordinate3D);
    Point3D getPoint3D(Point3D, Coordinate3D);
    Mesh1D createMesh1D(int, int *);
    Mesh3D createMesh3D(Point3D, int *, int *, int *);
    Mesh3D setMesh3D(char *);
    void setGrid(Mesh3D *, char *);
    int getNumberOfFine(int , int *);
    Mesh3D generateFineMesh(Mesh3D);
    Coordinate1D createCoordinate(Mesh1D, int, int);
    Coordinate3D createCoordinate3D(Mesh3D, Point3DDouble, Point3DDouble);
    Coordinate3D setCoordinate(char *);
    Coordinate3D change2Sphere(Coordinate3D, int);
    Point3D searchFineBase(Point3D, Coordinate3D);
    void copyMesh3D(Mesh3D *, Mesh3D *);
    """

    return header + func
