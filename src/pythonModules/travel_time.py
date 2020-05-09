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
        location = station.stationField.location
        tableField = _FDtomoC.lib.sphfd_exec(model.modelField, location)
        travel_time_table = TravelTimeTable(tableField = tableField)
        return travel_time_table
