#include "common/grid.h"
int getNumberOfXfine(Mesh mesh){
    int sum = 1;
    for(int i = 1; i < mesh.numberOfx; i++){
        sum += mesh.igridx[i - 1];
    }
    return sum;
}

int getNumberOfYfine(Mesh mesh){
    int sum = 1;
    for(int i = 1; i < mesh.numberOfy; i++){
        sum += mesh.igridy[i - 1];
    }
    return sum;
}

int getNumberOfZfine(Mesh mesh){
    int sum = 1;
    for(int i = 1; i < mesh.numberOfz; i++){
        sum += mesh.igridz[i - 1];
    }
    return sum;
}

int sizeofFine(Mesh mesh){
    return getNumberOfXfine(mesh) * getNumberOfYfine(mesh) * getNumberOfZfine(mesh);
}

float *getXFineMesh(Mesh mesh){
    int numberOfPoints = getNumberOfXfine(mesh);
    float *points = (float *)malloc(sizeof(float) * numberOfPoints);
    points[0] = mesh.origin.x;
    for(int i = 1; i < numberOfPoints; i++){
        points[i] = points[i - 1] + i * mesh.xSpace;
    }
    return points;
}
float *getYFineMesh(Mesh mesh){
    int numberOfPoints = getNumberOfYfine(mesh);
    float *points = (float *)malloc(sizeof(float) * numberOfPoints);
    points[0] = mesh.origin.y;
    for(int i = 1; i < numberOfPoints; i++){
        points[i] = points[i - 1] + i * mesh.ySpace;
    }
    return points;
}

float *getZFineMesh(Mesh mesh){
    int numberOfPoints = getNumberOfZfine(mesh);
    float *points = (float *)malloc(sizeof(float) * numberOfPoints);
    points[0] = mesh.origin.z;
    for(int i = 1; i < numberOfPoints; i++){
        points[i] = points[i - 1] + i * mesh.zSpace;
    }
    return points;
}

Point3D getFinePoint(Point3D point, Mesh mesh){
    float *gx = getXFineMesh(mesh);
    float *gy = getYFineMesh(mesh);
    float *gz = getZFineMesh(mesh);

    int x = (int)point.x;
    int y = (int)point.y;
    int z = (int)point.z;

    Point3D finePoint = {gx[x], gy[y], gz[z]};
    return finePoint;
}

float *getXMesh(Mesh mesh){
    float *points = (float *)malloc(sizeof(float) * mesh.numberOfx);
    points[0] = mesh.origin.x;
    for(int i = 1; i < mesh.numberOfx; i++){
        points[i] = points[i - 1] + mesh.igridx[i - 1] * mesh.xSpace;
    }
    return points;
}
float *getYMesh(Mesh mesh){
    float *points = (float *)malloc(sizeof(float) * mesh.numberOfy);
    points[0] = mesh.origin.y;
    for(int i = 1; i < mesh.numberOfy; i++){
        points[i] = points[i - 1] + mesh.igridy[i - 1] * mesh.ySpace;
    }
    return points;
}

float *getZMesh(Mesh mesh){
    float *points = (float *)malloc(sizeof(float) * mesh.numberOfz);
    points[0] = mesh.origin.z;
    for(int i = 1; i < mesh.numberOfz; i++){
        points[i] = points[i - 1] + mesh.igridz[i - 1] * mesh.zSpace;
    }
    return points;
}

Point3D getCoarsePoint(Point3D point, Mesh mesh){
    float *gx = getXMesh(mesh);
    float *gy = getYMesh(mesh);
    float *gz = getZMesh(mesh);

    int x = (int)point.x;
    int y = (int)point.y;
    int z = (int)point.z;

    Point3D CoarsePoint = {gx[x], gy[y], gz[z]};
    return CoarsePoint;
}

Mesh createMesh(SPEC spec){
    Mesh mesh;
	mesh.origin.x = spec.grid.x00;
    mesh.origin.y = spec.grid.y00;
    mesh.origin.z = spec.grid.z0;
    mesh.numberOfx = spec.grid.nxc;
    mesh.numberOfy = spec.grid.nyc;
    mesh.numberOfz = spec.grid.nzc;
    mesh.xSpace = mesh.ySpace = mesh.zSpace = spec.grid.h;

	mesh.igridx = (float *)malloc(sizeof(float)*spec.grid.nxc);
	mesh.igridy = (float *)malloc(sizeof(float)*spec.grid.nyc);
    mesh.igridz = (float *)malloc(sizeof(float)*spec.grid.nzc);
    memcpy(mesh.igridx, spec.grid.igridx, spec.grid.nxc * sizeof(float));
    memcpy(mesh.igridy, spec.grid.igridy, spec.grid.nyc * sizeof(float));
    memcpy(mesh.igridz, spec.grid.igridz, spec.grid.nzc * sizeof(float));

    return mesh;
}

void change2Sphere(Mesh mesh, int isElevation){
    double z0r;
    double rearth = 6371.0f, degrad = 0.017453292f, hpi = 1.570796f;
    
    mesh.origin.y *= degrad;
    if(isElevation){
        mesh.origin.y = hpi - glath(mesh.origin.y, mesh.origin.z, &z0r);
    } else {
        mesh.origin.y = hpi - glat(mesh.origin.y);
    }
    mesh.origin.x = mesh.origin.x * degrad;
    int space = mesh.xSpace; 
    mesh.ySpace = space / rearth / degrad;
    mesh.xSpace = fabs(space / (rearth * sin(mesh.origin.y)));
    mesh.xSpace = mesh.xSpace / degrad;

}

Point3D searchFineBase(Point3D point, Mesh mesh){
    float *gx = getXMesh(mesh);
    float *gy = getYMesh(mesh);
    float *gz = getZMesh(mesh);

    int i;
	for (i = 1; i < mesh.numberOfx; i++) {
		if (gx[i] > point.x)
			break;
	}
	if (i == mesh.numberOfx)
		i--;
	i--;

    int j;
	for (j = 1; j < mesh.numberOfy; j++) {
		if (gy[j] > point.y)
			break;
	}
	if (j == mesh.numberOfy)
		j--;
	j--;

    int k;
	for (k = 1; k < mesh.numberOfz; k++) {
		if (gz[k] > point.z)
			break;
	}
	if (k == mesh.numberOfz)
		k--;
	k--;

    Point3D base = {i, j, k};
    return base;
}
