import mesh, coordinate, velocity_model, environment, station, travel_time
import _FDtomoC

import abc

file_path = "../../data/small/FDtomo01.spec"
model1D_path = "../../data/small/TW_m30_mdl"
stafile = "../../data/small/runs_files/stationloc_out.txt";


loc_env = environment.LocEnv().create(file = file_path)
sphrayderv_env = environment.SphraydervEnv().create(file = file_path)
runlsqr_env = environment.RunlsqrEnv().create(file = file_path)

mesh3D = mesh.Mesh3D().create(file = file_path)
coarseCoordinate3D = coordinate.Coordinate3D().create(file = file_path)
fineCoordinate3D = coordinate.Coordinate3D().create(file = file_path)
fineMesh3D = mesh3D.generateFineMesh()
fineCoordinate3D.setMesh(fineMesh3D)

vpModel = velocity_model.VelocityModel1D()
vsModel = velocity_model.VelocityModel1D()
velocity_model.VelocityModel1D().setVelocityModel(model1D_path, vpModel, vsModel)

#tmp = [2 for i in range(21)]
#tmp[1:12] = [1 for i in range(11)]
tmp = [2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2]
fineMesh1D = mesh.Mesh1D().create(numberOfNode = 22, igrid = tmp)
fineCoordinate1D = coordinate.Coordinate1D().create(fineMesh1D, 2, -4)
fineVpModel = vpModel.transform(fineCoordinate1D)
CoarseModel3D = velocity_model.VelocityModel3D().create(coarseCoordinate3D, fineVpModel)
fineModel3D = CoarseModel3D.transform(fineCoordinate3D)

station_array = station.Station().createArray(file = stafile)
table = travel_time.TravelTimeTable().create(fineModel3D, station_array[2])