import abc
from mesh import Mesh3D
from point import Point
from cffi import FFI
import _FDtomoC

class TravelTimeTable(object):
    def __init__(self, name = None, mesh = None, time = None, tableField = None):
        self.name = name
        self.mesh = mesh
        self.time = time
        self.tableField = tableField

    def create(self, model, station):
        stationField = station.stationField
        travel_time_table = TravelTimeTable()
        travel_time_table.mesh = Mesh3D()
        travel_time_table.tableField = _FDtomoC.lib.sphfd_exec(model.modelField, stationField)
        travel_time_table.mesh.meshField = travel_time_table.tableField.mesh
        travel_time_table.mesh.numberOfNode.pointField = travel_time_table.tableField.mesh.numberOfNode
        travel_time_table.getClass()
        return travel_time_table

    def output(self, filename):
        tmp = _FDtomoC.ffi.new("char[]", filename)
        _FDtomoC.lib.outputTravelTime(self.tableField, tmp)

    def getClass(self):
        self.mesh.getClass()
        size = int(self.mesh.numberOfNode.x * self.mesh.numberOfNode.y * self.mesh.numberOfNode.z)
        self.time = _FDtomoC.ffi.unpack(self.tableField.time, size)
        self.name = _FDtomoC.ffi.string(self.tableField.name)

    def getField(self):
        meshField = self.mesh.getField()
        name = self.name
        time = _FDtomoC.ffi.new("float[]", self.time)
        #tableFieldPtr = _FDtomoC.ffi.new("travelTimeTable *", {'mesh' : meshField, 'name' : name, 'time' : time})
        tableField = _FDtomoC.lib.createTable(meshField, name, time)
        return tableField

    def removeField(self):
        self.tableField = None
        self.mesh.meshField = None
        self.mesh.numberOfNode.pointField = None
