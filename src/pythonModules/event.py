import abc
from cffi import FFI
import _FDtomoC
from derivative import Derivative
from residual import ResidualVector
from perturbation import Perturbation
class Event(object):
    def __init__(self, observation = None, earthquake = None, eventField = None):
        self.observation = observation
        self.earthquake = earthquake
        self.eventField = eventField

    def createArray(self, file = None):
        filename = _FDtomoC.ffi.new("char[]", file.encode('ascii'))
        eventField_list = _FDtomoC.lib.createEventList(filename)
        eventField_array = _FDtomoC.lib.EventList2Arr(eventField_list)
        
        listSize = _FDtomoC.lib.getEventCount(eventField_list)
        tmp = _FDtomoC.ffi.unpack(eventField_array, listSize)

        event_array = [Event(eventField = i) for i in tmp]

        return event_array

    def sphlocating(self, coordinate, table_array, locEnv, file):
        corField = coordinate.coordinateField
        tableFieldArray = [table_array[i].tableField for i in range(len(table_array))]
        tmp = _FDtomoC.ffi.new("travelTimeTable[]", tableFieldArray)
        locEnvField = locEnv.locEnvField

        filename = _FDtomoC.ffi.new("char[]", file.encode('ascii'))
        eventField_list = _FDtomoC.lib.createEventList(filename)

        table_size = len(tmp)
        eventFieldArray = _FDtomoC.lib.sphfdloc(corField, tmp, eventField_list, table_size, locEnvField)
        listSize = _FDtomoC.lib.getEventCount(eventField_list)
        event_array_tmp = _FDtomoC.ffi.unpack(eventFieldArray, listSize)
        event_array = [Event(eventField = i) for i in event_array_tmp]

        return event_array


    def singleLoc(self, coordinate, table_array, event, locEnv):
        corField = coordinate.coordinateField
        tableFieldArray = [table_array[i].tableField for i in range(len(table_array))]
        tmp = _FDtomoC.ffi.new("travelTimeTable[]", tableFieldArray)
        eventField = event.eventField
        locEnvField = locEnv.locEnvField

        newEventFielD = _FDtomoC.lib.singleLoc(corField, tmp, eventField, len(tmp), locEnvField)
        newEvent = Event(eventField = newEventFielD)

        return newEvent

    def sphRaytracing(self, model, table_array, event_array, event_size, station_array, table_size, sphrayderv_env):
        modelField = model.modelField
        tableFieldArray = [table_array[i].tableField for i in range(len(table_array))]
        table_tmp = _FDtomoC.ffi.new("travelTimeTable[]", tableFieldArray)

        eventFieldArray = [event_array[i].eventField for i in range(len(event_array))]
        event_tmp = _FDtomoC.ffi.new("Event[]", eventFieldArray)

        stationFieldArray = [station_array[i].stationField for i in range(len(station_array))]
        station_tmp = _FDtomoC.ffi.new("Station[]", stationFieldArray)

        sphraydervEnvField = sphrayderv_env.sphraydervEnvField
        commonEnvField = sphrayderv_env.commonEnvField

        data = _FDtomoC.lib.sphrayderv(modelField, table_tmp, event_tmp, event_size, station_tmp, table_size, sphraydervEnvField, commonEnvField)

        der = Derivative(data[0].mat)
        residual_vector = ResidualVector(data[0].b)
        return der, residual_vector

    def runlsqr(self, derivative, residual_vector, runlsqr_env):
        derivativeField = derivative.derivativeField
        residualField = residual_vector.residualField

        runlsqrEnvField = runlsqr_env.runlsqrEnvField
        commonEnvField = runlsqr_env.commonEnvField

        data = _FDtomoC.ffi.new("SPHRAYDERV_DATA *", {"mat" : derivativeField, "b" : residualField})
        perturbationField = _FDtomoC.lib.runlsqr(data, runlsqrEnvField, commonEnvField)
        perturbation = Perturbation(perturbationField)
        return perturbation
