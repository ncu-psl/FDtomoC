
import abc
from cffi import FFI
import _FDtomoC

class Coordinate(object):
    def __init__(self, mesh = None, origin = None,  space = None):
        self.mesh = mesh 
        self.space = space
        self.origin = origin
        self.coordinateField = None

    @abc.abstractmethod
    def create(self):
        return NotImplemented


class Coordinate1D(Coordinate):
    def create(self, mesh = None, space = None, origin = None, file = None):
        coordinate = Coordinate1D(mesh, space, origin)
        coordinate.coordinateField = _FDtomoC.lib.createCoordinate(mesh.meshField, space, origin)
        return  coordinate


class Coordinate3D(Coordinate):
    def create(self, mesh = None, origin = None, space = None, file = None):
        if (file != None):
            coordinate = Coordinate3D()
            tmp = _FDtomoC.ffi.new("char[]", file)
            coordinate.coordinateField = _FDtomoC.lib.setCoordinate(tmp)
            return coordinate

        coordinate = Coordinate3D(mesh, origin, space)
        coordinate.getField()
        return coordinate

    def setMesh(self, mesh):
        origin = self.coordinateField.origin
        space = self.coordinateField.space
        mesh = _FDtomoC.ffi.new("Mesh3D *", {'numberOfNode' : mesh.meshField.numberOfNode, \
                                              'gridx' : mesh.meshField.gridx, \
                                              'gridy' : mesh.meshField.gridy, \
                                              'gridz' : mesh.meshField.gridz  \
                                            })
        new_coordinateField = _FDtomoC.ffi.new("Coordinate3D *", {'mesh'   : mesh[0],   \
                                                                  'origin' : origin,    \
                                                                  'space'  : space      \
                                                                 })
        self.coordinateField = new_coordinateField[0]
        
    def getField(self):
        self.coordinateField = _FDtomoC.ffi.new("Coordinate3D *")
        self.coordinateField.origin = self.origin
        self.coordinateField.space = self.space
        self.setMesh(self.mesh)