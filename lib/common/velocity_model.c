#include "common/velocity_model.h"

void readVelocityModel1D(SPEC spec, velocityModel1D *vpModel, velocityModel1D *vsModel, char *interp){
	FILE *fp_one;
	fp_one = fopen(spec.onedfil, "r");
	if(!fp_one) {
        printf("Error on opening fp_one(%s)\n", spec.onedfil);
        assert(0);
	}
	vpModel->velocity = (float *)calloc(MAX1D, sizeof(float));
	vsModel->velocity = (float *)calloc(MAX1D, sizeof(float));
	char aline[MAXSTRLEN + 1];
	char pval[MAXSTRLEN + 1];

	double h, tmp = 0;
    int len, ierr;
	int ib = 0, ie = 0, lenv = 0, nvl = 0;
	int count = 0;
	get_line(fp_one, aline, &ierr);

	while(1){
	//---skip over header
	get_line(fp_one, aline, &ierr);
	if (ierr == 1)
		break;
	if (ierr != 0)
		continue;	
	ib = 0;
	get_field(fp_one, aline, ib, &ie, pval, &nvl, &ierr);
	sscanf(pval, "%lf", &h);
	ib = ie;
	get_field(fp_one, aline, ib, &ie, pval, &nvl, &ierr);
	float p = 0;
	sscanf(pval, "%f", &p);
	ib = ie;
	get_field(fp_one, aline, ib, &ie, pval, &nvl, &ierr);
	float s = 0;
	sscanf(pval, "%f", &s);
	ib = ie;
	get_field(fp_one, aline, ib, &ie, pval, &nvl, &ierr);

	if (count >= MAX1D) {
		printf(" Too many layers; maximum now is %d", MAX1D);
		assert(0);
	}
	vpModel->velocity[count] = p;
	if (spec.vs1d == 1) {
		vsModel->velocity[count] = s;
	} else {
		vsModel->velocity[count] = p/s;
	}
	
	if(count == 0){
		vpModel->coordinate.origin = h;
	}else{
		vec_push(&vpModel->coordinate.mesh1d.igrid, h - tmp);
	}
	tmp = h;
	interp[count] = pval[0];
	count++;
	}
	vpModel->coordinate.unit = 1;
	vpModel->coordinate.mesh1d.numberOfNode = count;
	vsModel->coordinate.mesh1d.numberOfNode = count;
	fclose(fp_one);
	
}

void transform1D(Coordinate1D coordinate, velocityModel1D *model, char *mode){
	int size = coordinate.mesh1d.numberOfNode;
	int vsize = model->coordinate.mesh1d.numberOfNode;
	float *points = getAxis(coordinate);
	float *vpoints = getAxis(model->coordinate);
	float *velocity = (float *)calloc(size, sizeof(float));
	model->coordinate = coordinate;
	velocity = linear_interpolation_array(points, vpoints, model->velocity, size, vsize, mode);
	free(model->velocity);
	model->velocity = velocity;
	free(points);
	free(vpoints);
}

velocityModel3D create3DModel(Coordinate3D coordinate, velocityModel1D model) {
	int meshSize3D = sizeOfMesh3D(coordinate.mesh);
	int axisSize1D = model.coordinate.mesh1d.numberOfNode;
	int axisSize3D = coordinate.mesh.numberOfNode.z;
	
	if(axisSize1D != axisSize3D){
		printf("Size Of coordinate mismatch with 1D velocity model!\n");
		assert(0);
	}
	/*
	//---unflatten the depths if required
    
	if (spec.iflat == 1) {
		for (int i = 0; i < spec.grid.nzc; i++) {
			gz[i] = uflatz(gz[i]);
		}
	}*/
//----generate the mode
	int index = 0;
	velocityModel3D model3D;
	int xsize = coordinate.mesh.numberOfNode.x;
	int ysize = coordinate.mesh.numberOfNode.y;
	int zsize = coordinate.mesh.numberOfNode.z;
	model3D.velocity = (float *)calloc(meshSize3D, sizeof(float));
	
	for (int k = 0; k < zsize; k++){
		for(int j = 0; j < ysize; j++){
			for(int i = 0; i < zsize; i++){
				model3D.velocity[index] = model.velocity[k];
				index++;
			}
		}
	}
	model3D.coordinate = coordinate;
	return model3D;
}

float getPointVel(Point3D point, velocityModel3D *model){
    int index = point.x * model->coordinate.mesh.numberOfNode.y * model->coordinate.mesh.numberOfNode.z +
                point.y * model->coordinate.mesh.numberOfNode.z +
                point.z;
    float vel = model->velocity[index];
    return vel;
}

Point3D getPoint3DModel(Point3D point, velocityModel3D *model){
	float *gx = getXAxis(model->coordinate);
    float *gy = getYAxis(model->coordinate);
    float *gz = getZAxis(model->coordinate);

    int x = (int)point.x;
    int y = (int)point.y;
    int z = (int)point.z;

    Point3D modelPoint = {gx[x], gy[y], gz[z]};
	free(gx);
    free(gy);
    free(gz);
    return modelPoint;
}

float trilinear_interpolation_base(Point3D point, Point3D base, velocityModel3D *model){
	Cell cells[2][2][2];
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			for(int k = 0; k < 2; k++){
				Point3D tmp = base;
				tmp.x += i;
				tmp.y += j;
				tmp.z += k;
				cells[i][j][k].point = getPoint3DModel(tmp, model);
				cells[i][j][k].value = getPointVel(tmp, model);
			}
		}
	}
	point = getPoint3D(point, model->coordinate);
	float vel = trilinear_interpolation(point, cells);
	return vel;
}



void transform3D(Coordinate3D coordinate, velocityModel3D *model){
	int meshSize3D = sizeOfMesh3D(coordinate.mesh);
	int vMeshSize3D = sizeOfMesh3D(model->coordinate.mesh);
	float xSize = coordinate.mesh.numberOfNode.x;
    float ySize = coordinate.mesh.numberOfNode.y;
    float zSize = coordinate.mesh.numberOfNode.z;
	float *velocity = (float *)calloc(meshSize3D, sizeof(float));

    int index = 0;
    for(int i = 0; i < xSize; i++){
        for(int j = 0; j < ySize; j++){
            for(int k = 0; k < zSize;k++){
                Point3D point = {i, j, k};
                Point3D location = getPoint3D(point, coordinate);
                Point3D base = searchFineBase(location, model->coordinate);
                float vel = trilinear_interpolation_base(point, base, model);
				velocity[index] = vel;
                index++;
            } 
        }
    }
	model->coordinate = coordinate;
	free(model->velocity);
	model->velocity = velocity;
	velocity = NULL;
}

velocityModel3D change2ColumnMajor(velocityModel3D model){
	velocityModel3D new_model;
	new_model.coordinate = model.coordinate;
	int x = model.coordinate.mesh.numberOfNode.x;
	int y = model.coordinate.mesh.numberOfNode.y;
	int z = model.coordinate.mesh.numberOfNode.z;
	float *vel = (float *)malloc(sizeof(float) * (x * y * z));
	int index = 0;
	for(int k = 0; k < z; k++){
		for(int j = 0; j < y; j++){
			for(int i = 0; i < x; i++){
				int offset = i * y * z +
							 j * y + 
							 k;
				vel[index++]=model.velocity[offset];
			}
		}
	}
	new_model.velocity = vel;
	return new_model;
}

void output3DModel(velocityModel3D model, char *filename){
	FILE *fp_tmp = fopen(filename, "wb");
	if(!fp_tmp){
		printf("Error happens while output model !\n");
		assert(0);
	}
	int size = sizeOfMesh3D(model.coordinate.mesh);
	fwrite(model.velocity, sizeof(float), size, fp_tmp);
	fclose(fp_tmp);
}

