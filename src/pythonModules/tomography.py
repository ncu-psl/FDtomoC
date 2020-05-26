from event import Event
from coordinate import Coordinate1D, Coordinate3D
from mesh import Mesh1D, Mesh3D
from velocity_model import VelocityModel1D, VelocityModel3D
from travel_time import TravelTimeTable
from perturbation import Perturbation
from residual import ResidualVector
from derivative import Derivative
from environment import CommonEnv, SphraydervEnv, LocEnv, RunlsqrEnv, MakenewmodEnv
import _FDtomoC
from mpi4py import MPI

class TomographyBuilder(object):
    def __init__(self):
        self.event_list = []
        self.station_list = []
        self.event_builder = None
        self.velocity_model_builder = None
        self.environment = None

    def Event(self, event = None):
        if(event != None):
            if(type(event) == list ):
                self.event_list.extend(event)
            else:
                self.event_list.append(event)
            return self
                
        if (self.event_builder != None):
            self.event_list.append(self.event_builder.getValue())
            
        self.event_builder = EventBuilder(self)
        return self.event_builder
    
    def Station(self, station):
        if(station != None):
            if type(station == list):
                self.station_list.extend(station)
            else:
                self.station_list.append(station)
        return self
    
    def VelocityModel(self, velocity_model = None):
        self.velocity_model_builder = VelocityModelBuilder(self)
        return self.velocity_model_builder

    def Environment(self, env):
        self.environment = env
        return self
    
    def execute(self, mode = None):
        if (self.event_builder != None):
            self.event_list.append(self.event_builder.getValue())

        loc_env = self.environment['loc_env']
        sphrayderv_env = self.environment['sphrayderv_env']
        runlsqr_env = self.environment['runlsqr_env']
        makenewmod_env = self.environment['makenewmod_env']
        
        coarseMesh3D = self.velocity_model_builder.coordinateBuilder.coarse_mesh
        fineMesh3D = self.velocity_model_builder.coordinateBuilder.fine_mesh
        origin = self.velocity_model_builder.coordinateBuilder.origin
        space = self.velocity_model_builder.coordinateBuilder.space
        
        numberOfNode = int(coarseMesh3D.meshField.numberOfNode.z)
        igrid = _FDtomoC.ffi.unpack(coarseMesh3D.meshField.gridz, numberOfNode - 1)
        zSpace = int(space.z)
        zOrigin = int(origin.z)
        coarseMesh1D = Mesh1D().create(numberOfNode = numberOfNode, igrid = igrid)
        cooarseCoordinate1D = Coordinate1D().create(coarseMesh1D, zSpace, zOrigin)
        
        vpModel1D = self.velocity_model_builder.vp_model
        vsModel1D = self.velocity_model_builder.vs_model
        coarseVpModel1D = vpModel1D.transform(cooarseCoordinate1D)
        coarseVsModel1D = vsModel1D.transform(cooarseCoordinate1D)
        
        coarseCoordinate3D = Coordinate3D().create(coarseMesh3D, space, origin)
        CoarseVpModel3D = VelocityModel3D().create(coarseCoordinate3D, coarseVpModel1D)
        CoarseVsModel3D = VelocityModel3D().create(coarseCoordinate3D, coarseVsModel1D)

        fineCoordinate3D = Coordinate3D().create(fineMesh3D, space, origin)
        fineVpModel3D = CoarseVpModel3D.transform(fineCoordinate3D)
        fineVsModel3D = CoarseVsModel3D.transform(fineCoordinate3D)

        if(mode == "MPI"):
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
            table_array = []

            for i in range(len(self.station_list)):
                if(rank == (i % size)):
                    table = TravelTimeTable().create(fineVpModel3D, self.station_list[i])
                    table.removeField()
                    table_array.append(table)

            table_array = comm.gather(table_array, root = 0)
            table_array = comm.bcast(table_array, root = 0)

            new_table_array = []
            for k in range(len(self.station_list)):
                station_name = _FDtomoC.ffi.string(self.station_list[k].stationField.name)
                for i in range(len(table_array)):
                    for j in range(len(table_array[i])):
                        if table_array[i][j].name == station_name:
                            table_array[i][j].tableField = table_array[i][j].getField()
                            new_table_array.append(table_array[i][j])

            new_event_array = []
            for i in range(len(self.event_list)):
                if (rank == (i % size)):
                    new_event = Event().singleLoc(fineCoordinate3D, new_table_array, self.event_list[i], loc_env)
                    new_event.removeField()
                    new_event_array.append(new_event)

            new_event_array = comm.gather(new_event_array, root = 0)


            event_array = []
            if (rank == 0):
                for i in range(len(new_event_array)):
                    for j in range(len(new_event_array[i])):
                        new_event_array[i][j].eventField = new_event_array[i][j].getField()
                        event_array.append(new_event_array[i][j])
                    
                event_size = len(event_array)
                table_size = len(new_table_array)
                derv, residual_vector = Event().sphRaytracing(CoarseVpModel3D, new_table_array, event_array, event_size, self.station_list, table_size, sphrayderv_env)
                perturbation = Event().runlsqr(derv, residual_vector, runlsqr_env)
                VelocityModel3D().makeNewModel(coarseCoordinate3D, CoarseVpModel3D, CoarseVsModel3D, perturbation, table_size, makenewmod_env)
        
        elif(mode == "Omp"):
            stationFieldArray = [self.station_list[i].stationField for i in range(len(self.station_list))]
            stationFieldArrayPtr = _FDtomoC.ffi.new("Station[]", stationFieldArray)

            tableArrayField = _FDtomoC.lib.sphfd(fineVpModel3D.modelField, stationFieldArrayPtr, len(stationFieldArray))

            eventFieldArray = [self.event_list[i].eventField for i in range(len(self.event_list))]
            eventFieldArrayPtr = _FDtomoC.ffi.new("Event[]", eventFieldArray)

            table_size = len(stationFieldArray)
            event_size = len(eventFieldArray)
            eventFieldArray = _FDtomoC.lib.sphfdloc(fineCoordinate3D.coordinateField, tableArrayField, table_size, eventFieldArrayPtr, event_size, loc_env.locEnvField)

            data = _FDtomoC.lib.sphrayderv(CoarseVpModel3D.modelField, tableArrayField, eventFieldArray, event_size, stationFieldArrayPtr, table_size, sphrayderv_env.sphraydervEnvField, sphrayderv_env.commonEnvField)
            perturbation = _FDtomoC.lib.runlsqr(data, runlsqr_env.runlsqrEnvField, runlsqr_env.commonEnvField)

            vpModelFieldPtr = _FDtomoC.ffi.new("velocityModel3D *", CoarseVpModel3D.modelField)
            vsModelFieldPtr = _FDtomoC.ffi.new("velocityModel3D *", CoarseVsModel3D.modelField)

            _FDtomoC.lib.makenewmod(coarseCoordinate3D.coordinateField, vpModelFieldPtr, \
                        vsModelFieldPtr, perturbation, table_size, makenewmod_env.makenewmodEnvField, makenewmod_env.commonEnvField)

        elif(mode == "Sequential"):
            table_list = []
            for i in range(len(self.station_list)):
                table = TravelTimeTable().create(fineVpModel3D, self.station_list[i])
                table_list.append(table)

            for i in range(len(self.station_list)):
                table = TravelTimeTable().create(fineVsModel3D, self.station_list[i])
                #table_list.append(table)
            
            new_event_list = []
            for i in range(len(self.event_list)):
                new_event = Event().singleLoc(fineCoordinate3D, table_list, self.event_list[i], loc_env)
                new_event_list.append(new_event)

            event_size = len(new_event_list)
            table_size = int(len(table_list))
            
            derv, residual_vector = Event().sphRaytracing(CoarseVpModel3D, table_list, new_event_list, event_size, self.station_list, table_size, sphrayderv_env)
            perturbation = Event().runlsqr(derv, residual_vector, runlsqr_env)
            VelocityModel3D().makeNewModel(coarseCoordinate3D, CoarseVpModel3D, CoarseVsModel3D, perturbation, table_size, makenewmod_env)
        
        return             
        
class EventBuilder():    
    def __init__(self, tomography_builder):
        self.observation_list = []
        self.earthquake = None
        self.earthquake_builder = None
        self.observation_builder = None
        self.tomography_builder = tomography_builder
        
    def __getattr__(self, name):
        def _method_missing(*args, **kwargs):
            if(name == 'Event'):
                return self.tomography_builder.Event(args[0])
            elif(name == 'Station'):
                return self.tomography_builder.Station(args[0])
        return _method_missing
            
    def Earthquake(self, earthquake):
        if(earthquake != None):
            self.earthquake = earthquake
            return self
        
        self.earthquake_builder = EarthquakeBuilder(self)
        return self.earthquake_builder
    
    def Observation(self, observation = None):
        if(observation != None):
            self.observation_list.append(observation)
            return self
        
        if (self.observation_builder != None):
            self.observation_list.append(self.observation_builder.getValue())
            
        self.observation_builder = ObservationBuilder(self)
        return self.observation_builder

    def getValue(self):
        self.hypocenter_list.append(self.observation_list)
        return self.hypocenter_list
    
    def execute(self):
        if (self.hypocenter_builder != None):
            self.hypocenter_list.append(self.hypocenter_builder.getValue())
        
        if (self.observation_builder != None):
            self.observation_list.append(self.observation_builder.getValue())
        
        self.tomography_builder.execute()
        
class VelocityModelBuilder():
    def __init__(self, tomography_builder):
        self.vp_model = None
        self.vs_model = None
        self.coordinate = None
        self.coordinateBuilder = None
        self.tomography_builder = tomography_builder
        
    def __getattr__(self, name):
        def _method_missing(*args, **kwargs):
            if(name == 'execute'):
                return self.tomography_builder.execute(*args, **kwargs)

        return _method_missing
        
    def ReferenceModel(self, vp_model, vs_model):
        self.vp_model = vp_model
        self.vs_model = vs_model
        return self

    def Coordinate(self, coordinate = None):
        self.coordinateBuilder = CoordinateBuilder(self)
        return  self.coordinateBuilder
    
    def getValue(self):
        if (self.coordinateBuilder != None):
            self.coordinate = self.coordinateBuilder.getValue()
        
        return [self.coordinate, self.reference_model]

class EarthquakeBuilder():
    def __init__(self, event_builder):
        self.event_builder = event_builder
        self.location = None
        self.time = None
        
    def __getattr__(self, name):
        def _method_missing(*args, **kwargs):
            if(name == 'Station'):
                return self.event_builder.Station(args[0])
        return _method_missing
        
    def Location(self, location):
        self.location = location
        return self
    
    def Time(self, time):
        self.time = time

    def Earthquake(self, time):
        return NotImplemented
        
    def getValue(self):
        self.hypocenter_list.append([self.time, self.location])
    
class ObservationBuilder():    
    def __init__(self, event_builder):
        self.observation_list = []
        self.time = None
        self.station = None
        self.setting = None
        self.event_builder = event_builder
        
    def __getattr__(self, name):
        def _method_missing(*args, **kwargs):
            if(name == 'Observation'):
                return self.event_builder.Observation(args[0])
        return _method_missing
        
    def Observation(self, observation = None):
        return self.event_builder.Observation()
    
    def Station(self, station):
        self.station = station
        return self
        
    def Time(self, time):
        self.time = time
        return self
        
    def Setting(self, setting):
        self.setting = setting
        return self

    def getValue(self):
        self.observation_list.append([self.station, self.time])
        return self.observation_list


class CoordinateBuilder():
    def __init__(self, velocity_model_builder):
        self.velocity_model_builder = velocity_model_builder
        self.coarse_mesh = None
        self.fine_mesh = None
        self.origin = None
        self.space = None

    def __getattr__(self, name):
        def _method_missing(*args, **kwargs):
            if(name == 'ReferenceModel'):
                return self.velocity_model_builder.ReferenceModel(args[0], args[1])

        return _method_missing    

    def Mesh(self, coarse_mesh, fine_mesh):
        self.coarse_mesh = coarse_mesh
        self.fine_mesh = fine_mesh
        return self
    
    def Origin(self, origin):
        self.origin = origin
        return self
    
    def Space(self, space):
        self.space = space
        return self

    def getValue(self):
        self.cooridinate.append([self.mesh, self.origin, self.space])
        return self.cooridinate


    
