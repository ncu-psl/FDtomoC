
from cffi import FFI
import _FDtomoC

class Coordinate(object):
    def __init__(self, mesh = None, space = None, origin = None):
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
    def create(self, mesh = None, space = None, origin = None, file = None):
        if (file != None):
            coordinate = Coordinate3D()
            tmp = _FDtomoC.ffi.new("char[]", file)
            coordinate.coordinateField = _FDtomoC.lib.setCoordinate(tmp)
            return coordinate
        
