import abc
from cffi import FFI
import _FDtomoC
from derivative import Derivative
from residual import ResidualVector
from perturbation import Perturbation
from point import Point

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

        event_array = []
        for i in range(len(tmp)):
            event = Event()
            event.earthquake = Earthquake()
            event.earthquake.location = Point()
            event.earthquake.time = Time()
            event.observation = Observation()
            event.observation.setting = Setting()

            event.eventField = tmp[i]
            event.earthquake.earthquakeField = tmp[i].earthquake
            event.earthquake.location.pointField = tmp[i].earthquake.location
            event.earthquake.time.timeField = tmp[i].earthquake.time


            event.getClass()
            event_array.append(event)

        return event_array

    def getClass(self):
        size = _FDtomoC.lib.getTimeCount(self.eventField.observedTimeList)

        self.observation.station_name = []
        self.observation.time_list = []
        currentTimeNode = self.eventField.observedTimeList
        for i in range(size):
            tmp = self.eventField.station_name_list[i]
            station_name = _FDtomoC.ffi.string(tmp)
            self.observation.station_name.append(station_name)

            timeField = currentTimeNode[0].time
            time = Time()
            time.timeField = timeField
            time.getClass()
            self.observation.time_list.append(time)

            currentTimeNode = currentTimeNode[0].next

        self.observation.setting.event_id = _FDtomoC.ffi.string(self.eventField.evid)
        self.observation.setting.phase = _FDtomoC.ffi.string(self.eventField.phase)
        self.observation.setting.isgood = _FDtomoC.ffi.unpack(self.eventField.isgood, size)
        self.observation.setting.rwts = _FDtomoC.ffi.unpack(self.eventField.rwts, size)

        self.earthquake.getClass()

    def getField(self):
        earthquakeField = self.earthquake.getField()
        evid = _FDtomoC.ffi.new("char[]", self.observation.setting.event_id)
        phase = _FDtomoC.ffi.new("char[]", self.observation.setting.phase)
        isgood = _FDtomoC.ffi.new("int[]", self.observation.setting.isgood)
        rwts = _FDtomoC.ffi.new("float[]", self.observation.setting.rwts)
        timeFieldArray = [time.getField() for time in self.observation.time_list]

        time_tmp = _FDtomoC.ffi.new("Time[]", timeFieldArray)
        timeFieldList = _FDtomoC.lib.TimeList2Arr(time_tmp, len(timeFieldArray))

        nameList = []
        for name in self.observation.station_name:
            station_name = _FDtomoC.ffi.new("char[133]", name)
            nameList.append(station_name)

        nameListField = _FDtomoC.ffi.new("char[maxobs][133]", nameList)
        EventNodeFieldPtr = _FDtomoC.lib.createEventNode(earthquakeField, nameListField, timeFieldList, phase, \
                                            rwts, isgood, evid)
        
        return EventNodeFieldPtr[0].event
 

    def removeField(self):
        self.eventField = None
        self.earthquake.earthquakeField = None
        self.earthquake.time.timeField = None
        self.earthquake.location.pointField = None
        for time in self.observation.time_list:
            time.timeField = None

    def sphlocating(self, coordinate, table_array, locEnv, event_array):
        corField = coordinate.coordinateField
        tableFieldArray = [table_array[i].tableField for i in range(len(table_array))]
        tableFieldArrayPtr = _FDtomoC.ffi.new("travelTimeTable[]", tableFieldArray)

        eventFieldArray = [event_array[i].eventField for i in range(len(event_array))]
        eventFieldArrayPtr = _FDtomoC.ffi.new("Event[]", eventFieldArray)

        locEnvField = locEnv.locEnvField

        table_size = len(tableFieldArray)
        event_size = len(eventFieldArray)
        eventFieldArray = _FDtomoC.lib.sphfdloc(corField, tableFieldArrayPtr, table_size, eventFieldArrayPtr, event_size, locEnvField)
        event_array_tmp = _FDtomoC.ffi.unpack(eventFieldArray, event_size)
        event_array = [Event(eventField = i) for i in event_array_tmp]

        return event_array


    def singleLoc(self, coordinate, table_array, event, locEnv):
        corField = coordinate.coordinateField
        tableFieldArray = [table_array[i].tableField for i in range(len(table_array))]
        tmp = _FDtomoC.ffi.new("travelTimeTable[]", tableFieldArray)
        eventField = event.eventField
        locEnvField = locEnv.locEnvField
        table_size = int(len(tmp)/2)

        newEventField = _FDtomoC.lib.singleLoc(corField, tmp, eventField, table_size, locEnvField)
        if(newEventField.evid[0] == b'0'):
            return -1
        else:
            newEvent = Event()
            newEvent.earthquake = Earthquake()
            newEvent.earthquake.location = Point()
            newEvent.earthquake.time = Time()
            newEvent.observation = Observation()
            newEvent.observation.setting = Setting()

            newEvent.eventField = newEventField
            newEvent.earthquake.earthquakeField = newEventField.earthquake
            newEvent.earthquake.location.pointField = newEventField.earthquake.location
            newEvent.earthquake.time.timeField = newEventField.earthquake.time

            newEvent.getClass()
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

class Earthquake(object):
    def __init__(self, location = None, time = None):
        self.location = location
        self.time = time
        self.earthquakeField = None

    def getClass(self):
        self.location.getClass()
        self.time.getClass()

    def getField(self):
        locationField = self.location.getField()
        timeField = self.time.getField()
        earthquakeFieldPtr = _FDtomoC.ffi.new("Earthquake *", {'location' : locationField, 'time' : timeField})
        return earthquakeFieldPtr[0]

class Observation(object):
    def __init__(self, station_name = None, time_list = None, setting = None):
        self.station_name = station_name
        self.time_list = time_list
        self.setting = setting
        self.observationField = None

    def getClass(self):
        self.time.getClass()
        station_name = _FDtomoC.ffi.new("char[]", self.station_name)
        #self.setting.event_id = _FDtomoC.ffi.string(self.event)

class Setting(object):
    def __init__(self, event_id = None, phase = None, isgood = None, rwts = None):
        self.event_id = event_id
        self.phase = phase
        self.isgood = isgood
        self.rwts = rwts        

class Time(object):
    def __init__(self, year = None, day = None, hour = None, minute = None, second = None):
        self.year = year
        self.day = day
        self.hour = hour
        self.minute = minute
        self.second = second
        self.timeField = None
    
    def getField(self):
        timeFieldPtr = _FDtomoC.ffi.new("Time *", {'iyr' : self.year, 'jday' : self.day, 'ihr' : self.hour, \
                                                    'imn' : self.minute, 'sec' : self.second})
        return timeFieldPtr[0]


    def getClass(self):
        self.year = self.timeField.iyr
        self.day = self.timeField.jday
        self.hour = self.timeField.ihr
        self.minute = self.timeField.imn
        self.second = self.timeField.sec