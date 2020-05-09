#include "common/velocity_model.h"

void readVelocityModel1D(char *model1D_path, velocityModel1D *vpModel, velocityModel1D *vsModel, char *interp){
	int vs1d = 1;
	FILE *fp_one;
	fp_one = fopen(model1D_path, "r");
	if(!fp_one) {
        printf("Error on opening fp_one(%s)\n", model1D_path);
        assert(0);
	}
	vpModel->velocity = (float *)calloc(MAX1D, sizeof(float));
	vsModel->velocity = (float *)calloc(MAX1D, sizeof(float));
	vpModel->coordinate.mesh.igrid = (int *)calloc(MAX1D, sizeof(int));
	vsModel->coordinate.mesh.igrid = (int *)calloc(MAX1D, sizeof(int));

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
	if (vs1d == 1) {
		vsModel->velocity[count] = s;
	} else {
		vsModel->velocity[count] = p/s;
	}
	
	if(count == 0){
		vpModel->coordinate.origin = h;
	}else{
		vpModel->coordinate.mesh.igrid[count - 1] = h - tmp;
	}
	tmp = h;
	interp[count] = pval[0];
	count++;
	}
	vpModel->coordinate.space = 1;
	vpModel->coordinate.mesh.numberOfNode = count;
	vsModel->coordinate.mesh.numberOfNode = count;
	fclose(fp_one);
	
}

velocityModel1D transform1D(Coordinate1D coordinate, velocityModel1D model, char *mode){
	velocityModel1D new_model;
	int size = coordinate.mesh.numberOfNode;
	int vsize = model.coordinate.mesh.numberOfNode;
	float *points = getAxis(coordinate);
	float *vpoints = getAxis(model.coordinate);
	float *velocity = (float *)calloc(size, sizeof(float));
	new_model.coordinate = coordinate;
	velocity = linear_interpolation_array(points, vpoints, model.velocity, size, vsize, mode);
	new_model.velocity = velocity;
	free(points);
	free(vpoints);
	return new_model;
}

velocityModel3D create3DModel(Coordinate3D coordinate, velocityModel1D model) {
	int meshSize3D = sizeOfMesh3D(coordinate.mesh);
	int axisSize1D = model.coordinate.mesh.numberOfNode;
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
			for(int i = 0; i < xsize; i++){
				model3D.velocity[index++] = model.velocity[k];
			}
		}
	}
	if(index > meshSize3D){
		printf("Error happens while getting velocity!\n");
		assert(0);
	}
	model3D.coordinate = coordinate;
	return model3D;
}

float getPointVel(Point3D point, velocityModel3D *model){
    int index = point.x +
                point.y * model->coordinate.mesh.numberOfNode.x +
                point.z * model->coordinate.mesh.numberOfNode.x * model->coordinate.mesh.numberOfNode.y;
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

float trilinear_interpolation_base(Point3D point, Point3D base, Coordinate3D coordinate, velocityModel3D *model){
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
	point = getPoint3D(point, coordinate);
	float vel = trilinear_interpolation2(point, cells);
	return vel;
}



velocityModel3D transform3D(Coordinate3D coordinate, velocityModel3D model){
	velocityModel3D new_model;
	int meshSize3D = sizeOfMesh3D(coordinate.mesh);
	int modelSize3D = sizeOfMesh3D(model.coordinate.mesh);
	int vMeshSize3D = sizeOfMesh3D(model.coordinate.mesh);
	float xSize = coordinate.mesh.numberOfNode.x;
    float ySize = coordinate.mesh.numberOfNode.y;
    float zSize = coordinate.mesh.numberOfNode.z;
	float *velocity = (float *)calloc(meshSize3D, sizeof(float));

	for(int i = 0; i < modelSize3D; i++){
		model.velocity[i] = 1. / model.velocity[i];
	}


    int index = 0;
    for(int k = 0; k < zSize; k++){
        for(int j = 0; j < ySize; j++){
            for(int i = 0; i < xSize; i++){
                Point3D point = {i, j, k};
                Point3D location = getPoint3D(point, coordinate);
                Point3D base = searchFineBase(location, model.coordinate);
                float vel = trilinear_interpolation_base(point, base, coordinate, &model);
				velocity[index] = 1.f / vel;
                index++;
            } 
        }
    }
	new_model.coordinate = coordinate;
	new_model.velocity = velocity;
	return new_model;
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

