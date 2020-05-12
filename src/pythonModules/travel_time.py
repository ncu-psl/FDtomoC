import abc
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
        tableField = _FDtomoC.lib.sphfd_exec(model.modelField, stationField)
        travel_time_table = TravelTimeTable(tableField = tableField)
        return travel_time_table

    def output(self, filename):
        tmp = _FDtomoC.ffi.new("char[]", filename)
        _FDtomoC.lib.outputTravelTime(self.tableField, tmp)
