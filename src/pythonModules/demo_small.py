import coordinate, environment, travel_time
from velocity_model import VelocityModel1D, VelocityModel3D
from mesh import Mesh1D, Mesh3D
from station import Station
from event import Event
from point import PointDouble
from travel_time import TravelTimeTable
from tomography import *
from mpi4py import MPI
import _FDtomoC

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

origin = PointDouble(121.30, 24.30, -2.0)

TomographyBuilder() \
    .Environment(environment) \
    .Event(event) \
    .Station(station) \
    .VelocityModel() \
        .Coordinate() \
            .Mesh(coarseMesh3D, fineMesh3D) \
            .Origin(origin) \
            .Space(PointDouble(2,2,2)) \
        .ReferenceModel(vpModel, vsModel) \
    .execute(mode = 'MPI', count = 1) 
