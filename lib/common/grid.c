#include "common/grid.h"
int sizeOfMesh3D(Mesh3D mesh){
    return mesh.numberOfNode.x * mesh.numberOfNode.y * mesh.numberOfNode.z;
}

float *getAxis(Coordinate1D coordinate){
    int numberOfPoints = coordinate.mesh.numberOfNode;
    float *points = (float *)malloc(sizeof(float) * numberOfPoints);
    points[0] = coordinate.origin;
    for(int i = 1; i < numberOfPoints; i++){
        points[i] = points[i - 1] + coordinate.space * coordinate.mesh.igrid[i - 1];
    }
    return points;
}


float *getXAxis(Coordinate3D coordinate){
    int numberOfPoints = coordinate.mesh.numberOfNode.x;
    float *points = (float *)calloc(numberOfPoints, sizeof(float));

    points[0] = coordinate.origin.x;
    for(int i = 1; i < numberOfPoints; i++){
        points[i] = points[i - 1] + coordinate.space.x * coordinate.mesh.gridx[i - 1];
    }
    return points;
}
float *getYAxis(Coordinate3D coordinate){
    int numberOfPoints = coordinate.mesh.numberOfNode.y;
    float *points = (float *)calloc(numberOfPoints, sizeof(float));
    points[0] = coordinate.origin.y;
    for(int i = 1; i < numberOfPoints; i++){
        points[i] = points[i - 1] + coordinate.space.y * coordinate.mesh.gridy[i - 1];
    }
    return points;
}

float *getZAxis(Coordinate3D coordinate){
    int numberOfPoints = coordinate.mesh.numberOfNode.z;
    float *points = (float *)calloc(numberOfPoints, sizeof(float));
    points[0] = coordinate.origin.z;
    for(int i = 1; i < numberOfPoints; i++){
        points[i] = points[i - 1] + coordinate.space.z * coordinate.mesh.gridz[i - 1];
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
    free(gx);
    free(gy);
    free(gz);
    return finePoint;
}

Mesh1D createMesh1D(int numberOfNode, int *igrid){
    Mesh1D mesh;
    mesh.numberOfNode = numberOfNode;
	mesh.igrid = calloc(numberOfNode - 1, sizeof(int));
    memcpy(mesh.igrid, igrid, sizeof(int) * (numberOfNode - 1));
    return mesh;
}

Mesh3D createMesh3D(Point3D numberOfNode, int *gridx, int *gridy, int *gridz){
    Mesh3D mesh;
    mesh.numberOfNode = numberOfNode;
    mesh.gridx = malloc((numberOfNode.x - 1) * sizeof(int));
    mesh.gridy = malloc((numberOfNode.y - 1) * sizeof(int));
    mesh.gridz = malloc((numberOfNode.z - 1) * sizeof(int));
    memcpy(mesh.gridx, gridx, sizeof(int) * (numberOfNode.x - 1));
    memcpy(mesh.gridy, gridy, sizeof(int) * (numberOfNode.y - 1));
    memcpy(mesh.gridz, gridz, sizeof(int) * (numberOfNode.z - 1));
    return mesh;
}

void copyMesh1D(Mesh1D *dest_mesh, Mesh1D *src_mesh){
    int size = src_mesh->numberOfNode;
    dest_mesh->numberOfNode = src_mesh->numberOfNode;
    dest_mesh->igrid = malloc(sizeof(int) * (size - 1));
    memcpy(dest_mesh->igrid, src_mesh->igrid, sizeof(int) * (size - 1));
}


Mesh3D setMesh3D(char *spec_file){
    Mesh3D mesh;
    char *mvals[3] = { "nxc\0", "nyc\0", "nzc\0"};
	char pval[MAXSTRLEN + 1];
    int len, ierr;

    FILE *fp_spc;
    fp_spc = fopen(spec_file, "r");
	if (!fp_spc) {
		printf("(Error in read_spec.c)read fp_spc file error.\n");
		assert(0);
	}

    int i;
	for (i = 0; i < 3; i++) {
		get_vars(fp_spc, mvals[i], pval, &len, &ierr);
		if (ierr == 1) {
			read_error(mvals[i], "variable", fp_spc);
		}
		if (i == 0) {
			sscanf(pval, "%f", &mesh.numberOfNode.x);
		} else if (i == 1) {
			sscanf(pval, "%f", &mesh.numberOfNode.y);
		} else if (i == 2) {
			sscanf(pval, "%f", &mesh.numberOfNode.z);
		}
    }

    mesh.gridx = malloc(sizeof(int) * (mesh.numberOfNode.x - 1));
    mesh.gridy = malloc(sizeof(int) * (mesh.numberOfNode.y - 1));
    mesh.gridz = malloc(sizeof(int) * (mesh.numberOfNode.z - 1));
    setGrid(&mesh, spec_file);
    fclose(fp_spc);
    return mesh;
}

void copyMesh3D(Mesh3D *dest_mesh, Mesh3D *src_mesh){
    int xSize = (int)src_mesh->numberOfNode.x;
    int ySize = (int)src_mesh->numberOfNode.y;
    int zSize = (int)src_mesh->numberOfNode.z;

    *dest_mesh = *src_mesh;
    dest_mesh->gridx = malloc(sizeof(int) * (xSize - 1));
    dest_mesh->gridy = malloc(sizeof(int) * (ySize - 1));
    dest_mesh->gridz = malloc(sizeof(int) * (zSize - 1));

    memcpy(dest_mesh->gridx, src_mesh->gridx, sizeof(int) * (xSize - 1));
    memcpy(dest_mesh->gridy, src_mesh->gridy, sizeof(int) * (ySize - 1));
    memcpy(dest_mesh->gridz, src_mesh->gridz, sizeof(int) * (zSize - 1));
}



void setGrid(Mesh3D *mesh, char *spec_file){
    FILE *fp_spc;
    fp_spc = fopen(spec_file, "r");
	if (!fp_spc) {
		printf("(Error in read_spec.c)read fp_spc file error.\n");
		assert(0);
	}

	char aline[MAXSTRLEN + 1], varname[MAXSTRLEN + 1], parval[MAXSTRLEN + 1];
	int len, ierr;
	int ib = 0, ie = 0, lenv = 0, nvl = 0;
	a11: get_line(fp_spc, aline, &ierr);
	aline[MAXSTRLEN] = '\0';
	if (ierr == 1)
		goto a12;
	if (ierr != 0)
		goto a11;
	get_field(fp_spc, aline, ib, &ie, varname, &lenv, &ierr);
	if (strncmp(varname, "igridx", lenv) != 0)
		goto a11;
	ib = ie;
	get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
	sscanf(parval, "%d", &mesh->gridx[0]);
	int k;
	for (k = 1; k < mesh->numberOfNode.x - 1; k++) {
		ib = ie;
		get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
		sscanf(parval, "%d", &mesh->gridx[k]);
	}
	a12: rewind(fp_spc);
	a13: get_line(fp_spc, aline, &ierr);
	if (ierr == 1)
		goto a14;
	if (ierr != 0)
		goto a13;
	ib = 0;
	get_field(fp_spc, aline, ib, &ie, varname, &lenv, &ierr);
	if (strncmp(varname, "igridy", lenv) != 0)
		goto a13;
	ib = ie;
	get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
	sscanf(parval, "%d", &mesh->gridy[0]);
	for (k = 1; k < mesh->numberOfNode.y - 1; k++) {
		ib = ie;
		get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
		sscanf(parval, "%d", &mesh->gridy[k]);
	}
	a14: rewind(fp_spc);

	a15: get_line(fp_spc, aline, &ierr);
	if (ierr != 1) {
		if (ierr != 0)
			goto a15;
		ib = 0;
		get_field(fp_spc, aline, ib, &ie, varname, &lenv, &ierr);
		if (strncmp(varname, "igridz", lenv) != 0)
			goto a15;
		ib = ie;
		get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
		sscanf(parval, "%d", &mesh->gridz[0]);
		for (k = 1; k < mesh->numberOfNode.z - 1; k++) {
			ib = ie;
			get_field(fp_spc, aline, ib, &ie, parval, &nvl, &ierr);
			sscanf(parval, "%d", &mesh->gridz[k]);
		}
	}
	fclose(fp_spc);

}

int getNumberOfFine(int numberOfNode, int *igrid){
    int sum = 1;
    for(int i = 0; i < numberOfNode - 1; i++){
        sum += igrid[i];
    }
    return sum;
}

Mesh3D generateFineMesh(Mesh3D mesh){
    Mesh3D new_mesh;
    new_mesh.numberOfNode.x = getNumberOfFine(mesh.numberOfNode.x, mesh.gridx);
    new_mesh.numberOfNode.y = getNumberOfFine(mesh.numberOfNode.y, mesh.gridy);
    new_mesh.numberOfNode.z = getNumberOfFine(mesh.numberOfNode.z, mesh.gridz);
    new_mesh.gridx = malloc(sizeof(int) * (new_mesh.numberOfNode.x - 1));
    new_mesh.gridy = malloc(sizeof(int) * (new_mesh.numberOfNode.y - 1));
    new_mesh.gridz = malloc(sizeof(int) * (new_mesh.numberOfNode.z - 1));

    for(int i = 0; i < new_mesh.numberOfNode.x - 1; i++){
        new_mesh.gridx[i] = 1;
    }
    for(int i = 0; i < new_mesh.numberOfNode.y - 1; i++){
        new_mesh.gridy[i] = 1;
    }
    for(int i = 0; i < new_mesh.numberOfNode.z - 1; i++){
        new_mesh.gridz[i] = 1;
    }
    return new_mesh;
}

Coordinate3D createCoordinate3D(Mesh3D mesh, Point3DDouble space, Point3DDouble origin){
    Coordinate3D coordinate;
    copyMesh3D(&coordinate.mesh, &mesh);
    coordinate.space = space;
    coordinate.origin = origin;
    return coordinate;
}


Coordinate1D createCoordinate(Mesh1D mesh, int space, int origin){
    Coordinate1D coordinate;
    copyMesh1D(&coordinate.mesh, &mesh);
    coordinate.space = space;
    coordinate.origin = origin;
    return coordinate;
}

Coordinate3D setCoordinate(char *spec_file){
    Coordinate3D coordinate;
    coordinate.mesh = setMesh3D(spec_file);

    double h;
    char *mvals[4] = { "x0\0", "y0\0", "z0\0", "h\0" };
    double *varList[4] = {&coordinate.origin.x, &coordinate.origin.y, &coordinate.origin.z, &h};

	char pval[MAXSTRLEN + 1];
    int len, ierr;
	sscanf(spec_file, "%s", spec_file);
    FILE *fp_spc;

    fp_spc = fopen(spec_file, "r");
	if (!fp_spc) {
		printf("(Error in read_spec.c)read fp_spc file error.\n");
		assert(0);
	}

    for(int i = 0; i < 4; i++){
        get_vars(fp_spc, mvals[i], pval, &len, &ierr);
		if (ierr == 1) {
			read_error(mvals[i], "variable", fp_spc);
		}            
        sscanf(pval, "%lf", varList[i]);
    }
    coordinate.space.x = coordinate.space.y = coordinate.space.z = h;
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
    coordinate.space.y = space / rearth;
    coordinate.space.x = fabs(space / (rearth * sin(coordinate.origin.y)));
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
    free(gx);
    free(gy);
    free(gz);
    return base;
}
