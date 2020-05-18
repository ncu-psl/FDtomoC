#!/usr/bin/env python
from cffi import FFI

def model_header():

    typedef = '''
    #define MAX1D 1000
    typedef struct velocityModel1D_{
        Coordinate1D coordinate;
        float *velocity;
    }velocityModel1D;

    typedef struct velocityModel3D_{
        Coordinate3D coordinate;
        float *velocity;
    }velocityModel3D;

     typedef struct{
        int n;
        float *x;
        int *jndx;
        float *se;
    }RUNLSQR_DATA;

    typedef struct{
        struct vhead head;
        ...;
        float *vn;

    }MAKENEWMOD_DATA;
    '''
    func = '''
    void readVelocityModel1D(char *, velocityModel1D *, velocityModel1D *, char *);
    velocityModel1D transform1D(Coordinate1D, velocityModel1D, char *);
    velocityModel3D create3DModel(Coordinate3D, velocityModel1D);
    velocityModel3D transform3D(Coordinate3D, velocityModel3D);
    float getPointVel(Point3D, velocityModel3D *);
    Point3D getPoint3DModel(Point3D, velocityModel3D *);
    float trilinear_interpolation_base(Point3D , Point3D, Coordinate3D, velocityModel3D *);
    velocityModel3D change2ColumnMajor(velocityModel3D);
    void output3DModel(velocityModel3D, char *);
    void freeModel3D(velocityModel3D *);
    void freeModel1D(velocityModel1D *);
    velocityModel1D *createModel1D();
    void makenewmod(Coordinate3D, velocityModel3D *, velocityModel3D *, RUNLSQR_DATA *, int, MakenewmodEnv, CommonEnv);
    '''

    return typedef + func