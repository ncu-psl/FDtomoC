#include "common/grid.h"
int sizeOfMesh3D(Mesh3D mesh){
    return mesh.numberOfNode.x * mesh.numberOfNode.y * mesh.numberOfNode.z;
}

float *getXAxis(Coordinate3D coordinate){
    int numberOfPoints = coordinate.mesh.numberOfNode.x;
    float *points = (float *)malloc(sizeof(float) * numberOfPoints);
    points[0] = coordinate.origin.x;
    for(int i = 1; i < numberOfPoints; i++){
        points[i] = points[i - 1] + i * coordinate.mesh.space.x;
    }
    return points;
}
float *getYAxis(Coordinate3D coordinate){
    int numberOfPoints = coordinate.mesh.numberOfNode.y;
    float *points = (float *)malloc(sizeof(float) * numberOfPoints);
    points[0] = coordinate.origin.y;
    for(int i = 1; i < numberOfPoints; i++){
        points[i] = points[i - 1] + i * coordinate.mesh.space.y;
    }
    return points;
}

float *getZAxis(Coordinate3D coordinate){
    int numberOfPoints = coordinate.mesh.numberOfNode.z;
    float *points = (float *)malloc(sizeof(float) * numberOfPoints);
    points[0] = coordinate.origin.z;
    for(int i = 1; i < numberOfPoints; i++){
        points[i] = points[i - 1] + i * coordinate.mesh.space.z;
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
