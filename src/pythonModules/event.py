import abc
from cffi import FFI
import _FDtomoC

class Event(object):
    def __init__(self, observation = None, earthquake = None, eventField = None):
        self.observation = observation
        self.earthquake = earthquake
        self.eventField = eventField

    def createArray(self, file = None):
        filename = _FDtomoC.ffi.new("char[]", file)
        eventField_list = _FDtomoC.lib.createEventList(filename)
        eventField_array = _FDtomoC.lib.EventList2Arr(eventField_list)
        
        listSize = _FDtomoC.lib.getEventCount(eventField_list)
        tmp = _FDtomoC.ffi.unpack(eventField_array, listSize)

        event_array = [Event(eventField = i) for i in tmp]

        return event_array

    def singleLoc(self, coordinate, table_array, event, locEnv):
        corField = coordinate.coordinateField
        tableField_array = [table_array[i].tableField for i in range(len(table_array))]
        tmp = _FDtomoC.ffi.new("travelTimeTable[]", tableField_array)
        eventField = event.eventField
        locEnvField = locEnv.locEnvField

        newEventFielD = _FDtomoC.lib.singleLoc(corField, tmp, eventField, len(tmp), locEnvField)
        newEvent = Event(eventField = newEventFielD)

        return newEvent