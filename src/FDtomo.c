#include "common/read_spec.h"
#include "common/velocity_model.h"
#include "common/station.h"
#include "FDtomo/make1d.h"
#include "FDtomo/c2f.h"
#include "FDtomo/sphfd.h"
#include "FDtomo/sphfdloc.h"
#include "FDtomo/sphrayderv.h"
#include "FDtomo/runlsqr.h"
#include "FDtomo/makenewmod.h"
#include "common/grid.h"
#include "common/shared_variables.h"
#include "common/vec/vec.h"
#include <stdio.h>
int main(int argc, char *argv[]){
	
	char spec_path[MAXSTRLEN] = "../data/small/FDtomo01.spec";
	char *model1D_path = "../data/small/TW_m30_mdl";
	char *stafile = "../data/small/runs_files/stationloc_out.txt";
	char *leqsfil = "../data/small/runs_files/arrivals/All.txt";
	
	/*
	char spec_path[MAXSTRLEN] = "../data/HL_2km/run_spec_file01";
	char *model1D_path = "../data/HL_2km/CWB_mdl";
	char *stafile = "../data/HL_2km/run_files/HL171819.coords";
	char *leqsfil = "../data/HL_2km/run_files/arrivals/hyp.Fdtomo171819";
	*/

	CommonEnv common_env = setCommonEnv(spec_path);
	LocEnv loc_env = setLocEnv(spec_path);
	SphraydervEnv sphrayderv_env = setSphraydervEnv(spec_path);
	RunlsqrEnv runlsqr_env = setRunlsqrEnv(spec_path);
	MakenewmodEnv makenewmod_env = setMakeNewmodEnv(spec_path);
	Mesh3D mesh = setMesh3D(spec_path);
	Coordinate3D coarseCoordinate =  setCoordinate(spec_path);
	Coordinate3D fineCoordinate = setCoordinate(spec_path);

	Mesh3D fineMesh = generateFineMesh(mesh);
	copyMesh3D(&fineCoordinate.mesh, &fineMesh);
	
	velocityModel1D vpModel; 
	velocityModel1D vsModel; 
	char interp[MAX1D];
	readVelocityModel1D(model1D_path, &vpModel, &vsModel, interp);
	
	Mesh1D coarseMesh = createMesh1D(coarseCoordinate.mesh.numberOfNode.z, coarseCoordinate.mesh.gridz);
	Coordinate1D coarseCoordinate1D =  createCoordinate(coarseMesh, coarseCoordinate.space.z, coarseCoordinate.origin.z);
	
	velocityModel1D newVpModel = transform1D(coarseCoordinate1D, vpModel, interp);
	velocityModel3D coarseVpModel = create3DModel(coarseCoordinate, newVpModel);

	EventNode *event_list = createEventList(leqsfil);
	StationNode *station_list = createStationList(stafile, 1);
	int table_size = getStationCount(station_list);
	Station *station_array = StationList2Arr(station_list);
	Event *event_array = EventList2Arr(event_list);

	for(int i = 0; i < 1; i++){
		velocityModel3D fineVpModel = transform3D(fineCoordinate, coarseVpModel);
		
		velocityModel1D newVsModel = transform1D(coarseCoordinate1D, vsModel, interp);
		velocityModel3D coarseVsModel = create3DModel(coarseCoordinate, newVsModel);
		velocityModel3D fineVsModel = transform3D(fineCoordinate, coarseVsModel);

		travelTimeTable *vpTable_array = sphfdAll(fineVpModel, fineVsModel, station_array, table_size);
		
		int event_size = getEventCount(event_list);
		EventNode *new_event_list = sphfdloc(fineCoordinate, vpTable_array, table_size, event_array, event_size, loc_env);

		int new_event_size = getEventCount(new_event_list);
		Event *new_event_array = EventList2Arr(new_event_list);
		SPHRAYDERV_DATA *SPHRAYDERV = sphrayderv(coarseVpModel, vpTable_array, new_event_array, new_event_size, station_array, table_size, sphrayderv_env, common_env);

		RUNLSQR_DATA *RUNLSQR = runlsqr(SPHRAYDERV, runlsqr_env, common_env);
		makenewmod(coarseCoordinate, &coarseVpModel, &coarseVsModel, RUNLSQR, table_size, makenewmod_env, common_env);
		free(SPHRAYDERV);
		free(RUNLSQR);
		free(new_event_array);
		free(vpTable_array);
	}
	return 0;

}
