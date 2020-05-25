import abc
from cffi import FFI
from mesh import Mesh1D, Mesh3D
from point import PointDouble
import _FDtomoC

class Coordinate(object):
    def __init__(self, mesh = None, space = None,  origin = None):
        self.mesh = mesh 
        self.space = space
        self.origin = origin
        self.coordinateField = None

    @abc.abstractmethod
    def create(self):
        return NotImplemented

    @abc.abstractmethod
    def getField(self):
        return NotImplemented

    @abc.abstractmethod
    def getClass(self):
        return NotImplemented

class Coordinate1D(Coordinate):
    def create(self, mesh = None, space = None, origin = None, file = None):
        coordinate = Coordinate1D(mesh, space, origin)
        coordinate.coordinateField = coordinate.getField()
        coordinate.mesh.meshField = coordinate.coordinateField.mesh
        return  coordinate

    def getField(self):
        meshField = self.mesh.getField()
        coordinateField = _FDtomoC.lib.createCoordinate(meshField, self.space, self.origin)
        '''
        coordinateFieldPtr = _FDtomoC.ffi.new("Coordinate1D *", {'mesh'   : meshField, \
                                                                 'origin' : self.origin, \
                                                                 'space'  : self.space \
                                                                })
        '''
        return coordinateField


    def getClass(self):
        self.mesh.getClass()
        self.origin = self.coordinateField.origin
        self.space = self.coordinateField.space


class Coordinate3D(Coordinate):
    def __init__(self, mesh = None, space = None,  origin = None):
        self.mesh = mesh 
        self.space = space
        self.origin = origin
        self.coordinateField = None

    def create(self, mesh = None, space = None, origin = None, file = None):
        if (file != None):
            coordinate = Coordinate3D()
            coordinate.mesh = Mesh3D()
            coordinate.space = PointDouble()
            coordinate.origin = PointDouble()
            tmp = _FDtomoC.ffi.new("char[]", file.encode('ascii'))
            coordinate.coordinateField = _FDtomoC.lib.setCoordinate(tmp)
            coordinate.mesh.meshField = coordinate.coordinateField.mesh
            coordinate.mesh.numberOfNode.pointField = coordinate.coordinateField.mesh.numberOfNode
            coordinate.space.pointField = coordinate.coordinateField.space
            coordinate.origin.pointField = coordinate.coordinateField.origin
            coordinate.getClass()
            return coordinate

        coordinate = Coordinate3D(mesh, space, origin)
        coordinate.coordinateField = coordinate.getField()
        return coordinate

    def setMesh(self, mesh):
        origin = self.coordinateField.origin
        space = self.coordinateField.space
        #meshField = mesh.getField()  cause error that value of gridz cannot keep permanently
        meshFieldPtr = _FDtomoC.ffi.new("Mesh3D *", {'numberOfNode' : mesh.meshField.numberOfNode, \
                                        'gridx' : mesh.meshField.gridx, \
                                        'gridy' : mesh.meshField.gridy, \
                                        'gridz' : mesh.meshField.gridz  \
                                    })
        new_coordinateField = _FDtomoC.ffi.new("Coordinate3D *", {'mesh'   : meshFieldPtr[0], \
                                                                  'origin' : origin,    \
                                                                  'space'  : space      \
                                                                 })
        self.coordinateField = new_coordinateField[0]
        self.mesh.meshField = self.coordinateField.mesh
        self.space.pointField = self.coordinateField.space
        self.origin.pointField = self.coordinateField.origin
        
    def getField(self):
        meshField = self.mesh.getField()
        originField = self.origin.getField()
        spaceField = self.space.getField()
        coordinateField = _FDtomoC.lib.createCoordinate3D(meshField, spaceField, originField)
        '''
        coordinateFieldPtr = _FDtomoC.ffi.new("Coordinate3D *", {'mesh'   : meshField, \
                                                                 'origin' : originField, \
                                                                 'space'  : spaceField \
                                                                })
        '''
        return coordinateField

    def getClass(self):
        self.mesh.getClass()
        self.origin.getClass()
        self.space.getClass()