import abc
from cffi import FFI
import _FDtomoC

class Station(object):
    def __init__(self, name = None, location = None, stationField = None):
        self.name = name
        self.location = location
        self.stationField = stationField

    def create(self, name = None, location = None):
        station = Station(name, location)
        locationField = _FDtomoC.ffi.new("Point3DDouble *", location)
        nameField = _FDtomoC.ffi.new("char[]", name)
        stationField = _FDtomoC.ffi.new("Station *", {'name' : nameField, 'location' : locationField})
        station.stationField = stationField
        return station

    def createArray(self, station = None, file = None):
        if(file != None):
            filename = _FDtomoC.ffi.new("char[]", file.encode('ascii'))
            stationField_list = _FDtomoC.lib.createStationList(filename, 1)
            stationField_array = _FDtomoC.lib.StationList2Arr(stationField_list)

            stationSize = _FDtomoC.lib.getStationCount(stationField_list)
            tmp = _FDtomoC.ffi.unpack(stationField_array, stationSize)

            station_array = [Station(stationField = i) for i in tmp]
            
            return station_array

