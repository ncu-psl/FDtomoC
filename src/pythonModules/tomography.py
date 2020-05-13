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
    
    def execute(self):
        if (self.event_builder != None):
            self.event_list.append(self.event_builder.getValue())
        
        coarseMesh3D = self.velocity_model_builder.coordinateBuilder.coarse_mesh
        fineMesh3D = self.velocity_model_builder.coordinateBuilder.fine_mesh
        origin = self.velocity_model_builder.coordinateBuilder.origin
        space = self.velocity_model_builder.coordinateBuilder.space

        numberOfNode = int(coarseMesh3D.meshField.numberOfNode.z)
        igrid = _FDtomoC.ffi.unpack(coarseMesh3D.meshField.gridz, numberOfNode - 1)
        coarseMesh1D = Mesh1D().create(numberOfNode = numberOfNode, igrid = igrid)
        cooarseCoordinate1D = Coordinate1D().create(coarseMesh1D, 2, -4)

        vpModel1D = self.velocity_model_builder.vp_model
        vsModel1D = self.velocity_model_builder.vs_model
        coarseVpModel1D = vpModel1D.transform(cooarseCoordinate1D)
        coarseVsModel1D = vsModel1D.transform(cooarseCoordinate1D)

        coarseCoordinate3D = Coordinate3D().create(coarseMesh3D, origin, space)
        CoarseVpModel3D = VelocityModel3D().create(coarseCoordinate3D, coarseVpModel1D)
        CoarseVsModel3D = VelocityModel3D().create(coarseCoordinate3D, coarseVsModel1D)
        
        fineCoordinate3D = Coordinate3D().create(fineMesh3D, origin, space)
        fineVpModel3D = CoarseVpModel3D.transform(fineCoordinate3D)
        fineVsModel3D = CoarseVsModel3D.transform(fineCoordinate3D)

        table_list = []
        print(self.station_list[0].stationField.location.x)
        print(fineVpModel3D.modelField.coordinate.mesh.numberOfNode.x)
        for i in range(len(self.station_list)):
            table = TravelTimeTable().create(fineVpModel3D, self.station_list[i])
            table_list.append(table)


        loc_env = self.environment['loc_env']
        new_event_list = []
        for i in range(len(self.event_list)):
            new_event = Event().singleLoc(fineCoordinate3D, table_list, self.event_list[i], loc_env)
            new_event_list.append(new_event)


        sphrayderv_env = self.environment['sphrayderv_env']
        runlsqr_env = self.environment['runlsqr_env']
        makenewmod_env = self.environment['makenewmod_env']
        event_size = len(new_event_list)
        table_size = len(table_list)
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
                return self.tomography_builder.execute()

        return _method_missing
        
    def ReferenceModel(self, vp_model, vs_model):
        self.vp_model = vp_model
        self.vs_model = vs_model
        return self

    def Coordinate(self, coordinate = None):
        if (coordinate != None):
            self.coordinate = coordinate
            return self
        
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
        self.cooridinate = []
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


    
