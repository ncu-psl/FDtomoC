#include "common/grid.h"
int sizeOfMesh3D(Mesh3D mesh){
    return mesh.numberOfNode.x * mesh.numberOfNode.y * mesh.numberOfNode.z;
}

float *getXAxis(Coordinate3D coordinate){
    int numberOfPoints = coordinate.mesh.numberOfNode.x;
    float *points = (float *)malloc(sizeof(float) * numberOfPoints);
    points[0] = coordinate.origin.x;
    for(int i = 1; i < numberOfPoints; i++){
        points[i] = points[i - 1] + i * coordinate.mesh.xspace;
    }
    return points;
}
float *getYAxis(Coordinate3D coordinate){
    int numberOfPoints = coordinate.mesh.numberOfNode.y;
    float *points = (float *)malloc(sizeof(float) * numberOfPoints);
    points[0] = coordinate.origin.y;
    for(int i = 1; i < numberOfPoints; i++){
        points[i] = points[i - 1] + i * coordinate.mesh.yspace;
    }
    return points;
}

float *getZAxis(Coordinate3D coordinate){
    int numberOfPoints = coordinate.mesh.numberOfNode.z;
    float *points = (float *)malloc(sizeof(float) * numberOfPoints);
    points[0] = coordinate.origin.z;
    for(int i = 1; i < numberOfPoints; i++){
        points[i] = points[i - 1] + i * coordinate.mesh.zspace;
    }
    return points;
}

Point3D getPoint3D(Point3D point, Coordinate3D coordinate){
    float *gx = getXAxis(coordinate);
    float *gy = getYAxis(coordinate);
    float *gz = getZAxis(coordinate);

    int x = (int)point.x;
    int y = (int)point.y;
    int z = (int)point.z;

    Point3D finePoint = {gx[x], gy[y], gz[z]};
    return finePoint;
}

Mesh1D createMesh1D(int numberOfNode, int space, int *igrid){
    Mesh1D mesh;
    vec_init(&mesh.igrid);
    mesh.numberOfNode = numberOfNode;
    mesh.space = space;
    return mesh;
}

Mesh3D readFineMesh3D(SPEC spec){
    Mesh3D mesh;
    mesh.numberOfNode.x = getNumberOfXfine(spec);
    mesh.numberOfNode.y = getNumberOfYfine(spec);
    mesh.numberOfNode.z = getNumberOfZfine(spec);
    mesh.xspace = spec.grid.xSpace;
    mesh.yspace = spec.grid.ySpace;
    mesh.zspace = spec.grid.zSpace;
    vec_init(&mesh.igridx);
    vec_init(&mesh.igridy);
    vec_init(&mesh.igridz);

    for(int i = 1; i < spec.grid.nxc; i++){
        vec_push(&mesh.igridx, spec.grid.igridx[i - 1]);
    }
    for(int i = 1; i < spec.grid.nyc; i++){
        vec_push(&mesh.igridy, spec.grid.igridy[i - 1]);
    }
    for(int i = 1; i < spec.grid.nzc; i++){
        vec_push(&mesh.igridz, spec.grid.igridz[i - 1]);
    }

    return mesh;
}

Mesh3D readCoarseMesh3D(SPEC spec){
    Mesh3D mesh;
    mesh.numberOfNode.x = spec.grid.nxc;
    mesh.numberOfNode.y = spec.grid.nyc;
    mesh.numberOfNode.z = spec.grid.nzc;
}

Coordinate3D readFineCoordinate(SPEC spec){
    Coordinate3D coordinate;
    coordinate.mesh = readFineMesh3D(spec);
    coordinate.origin = (Point3DDouble){spec.grid.x00, spec.grid.y00, spec.grid.z0};
    coordinate.space = (Point3DDouble){spec.grid.xSpace, spec.grid.ySpace, spec.grid.zSpace};
    return coordinate;
}

Coordinate3D readCoarseCoordinate(SPEC spec){
    Coordinate3D coordinate;
    coordinate.mesh = readCoarseMesh3D(spec);
    coordinate.origin = (Point3DDouble){spec.grid.x00, spec.grid.y00, spec.grid.z0};
    return coordinate;
}

Coordinate3D change2Sphere(Coordinate3D coordinate, int isElevation){
    double z0r;
    double rearth = 6371.0f, degrad = 0.017453292f, hpi = 1.570796f;
    coordinate.origin.y *= degrad;
    if(isElevation){
        coordinate.origin.y = hpi - glath(coordinate.origin.y, coordinate.origin.z, &z0r);
    } else {
        coordinate.origin.y = hpi - glat(coordinate.origin.y);
    }
    coordinate.origin.x = coordinate.origin.x * degrad;
    int space = coordinate.space.x; 
    coordinate.space.y = space / rearth / degrad;
    coordinate.space.x = fabs(space / (rearth * sin(coordinate.origin.y)));
    coordinate.space.x = coordinate.space.x / degrad;
    return coordinate;
}

Point3D searchFineBase(Point3D point, Coordinate3D coordinate){
    float *gx = getXAxis(coordinate);
    float *gy = getYAxis(coordinate);
    float *gz = getZAxis(coordinate);

    int i;
	for (i = 1; i < coordinate.mesh.numberOfNode.x; i++) {
		if (gx[i] > point.x)
			break;
	}
	if (i == coordinate.mesh.numberOfNode.x)
		i--;
	i--;

    int j;
	for (j = 1; j < coordinate.mesh.numberOfNode.y; j++) {
		if (gy[j] > point.y)
			break;
	}
	if (j == coordinate.mesh.numberOfNode.y)
		j--;
	j--;

    int k;
	for (k = 1; k < coordinate.mesh.numberOfNode.z; k++) {
		if (gz[k] > point.z)
			break;
	}
	if (k == coordinate.mesh.numberOfNode.z)
		k--;
	k--;

    Point3D base = {i, j, k};
    return base;
}
