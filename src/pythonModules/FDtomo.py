import coordinate, velocity_model, environment, station, travel_time, event
from velocity_model import VelocityModel1D, VelocityModel3D
from mesh import Mesh1D, Mesh3D
from station import Station
from event import Event
import _FDtomoC
import abc
from tomography import *
file_path = "../../data/small/FDtomo01.spec"
model1D_path = "../../data/small/TW_m30_mdl"
stafile = "../../data/small/runs_files/stationloc_out.txt"
leqsfil = "../../data/small/runs_files/arrivals/All.txt"

loc_env = environment.LocEnv().create(file = file_path)
sphrayderv_env = environment.SphraydervEnv().create(file = file_path)
runlsqr_env = environment.RunlsqrEnv().create(file = file_path)
makenewmod_env = environment.MakenewmodEnv().create(file = file_path)

environment = {'loc_env' : loc_env, 'sphrayderv_env' : sphrayderv_env, 'runlsqr_env' : runlsqr_env, 'makenewmod_env' : makenewmod_env}

coarseMesh3D = Mesh3D().create(file = file_path)
fineMesh3D = coarseMesh3D.generateFineMesh()

event = Event().createArray(leqsfil)
station = Station().createArray(file = stafile)

vpModel = VelocityModel1D()
vsModel = VelocityModel1D()
VelocityModel1D().setVelocityModel(model1D_path, vpModel, vsModel)

TomographyBuilder() \
    .Environment(environment) \
    .Event(event) \
    .Station(station) \
    .VelocityModel() \
        .Coordinate() \
            .Mesh(coarseMesh3D, fineMesh3D) \
            .Origin([120.9, 23.8, -4.0]) \
            .Space([2,2,2]) \
        .ReferenceModel(vpModel, vsModel) \
    .execute() 

'''
loc_env = environment.LocEnv().create(file = file_path)
sphrayderv_env = environment.SphraydervEnv().create(file = file_path)
runlsqr_env = environment.RunlsqrEnv().create(file = file_path)
makenewmod_env = environment.MakenewmodEnv().create(file = file_path)


mesh3D = Mesh3D().create(file = file_path)
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
fineVsModel = vsModel.transform(fineCoordinate1D)

CoarseVpModel3D = velocity_model.VelocityModel3D().create(coarseCoordinate3D, fineVpModel)
fineVpModel3D = CoarseVpModel3D.transform(fineCoordinate3D)

CoarseVsModel3D = velocity_model.VelocityModel3D().create(coarseCoordinate3D, fineVsModel)


station_array = station.Station().createArray(file = stafile)
table_array = []
for i in range(len(station_array)):
   table = travel_time.TravelTimeTable().create(fineVpModel3D, station_array[i])
   table_array.append(table)

event_array = event.Event().createArray(leqsfil)
#for i in range(len(event_array)):
 #   new_event = event.Event().singleLoc(fineCoordinate3D, table_array, event_array[i], loc_env)
  #  new_event_array.append(new_event)

new_event_array = event.Event().sphlocating(fineCoordinate3D, table_array, loc_env, leqsfil)

event_size = len(new_event_array)
table_size = len(table_array)
derv, residual_vector = event.Event().sphRaytracing(CoarseVpModel3D, table_array, new_event_array, event_size, station_array, table_size, sphrayderv_env)
perturbation = event.Event().runlsqr(derv, residual_vector, runlsqr_env)
velocity_model.VelocityModel3D().makeNewModel(coarseCoordinate3D, CoarseVpModel3D, CoarseVsModel3D, perturbation, table_size, makenewmod_env)
'''