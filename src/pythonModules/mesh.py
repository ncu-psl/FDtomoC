import abc
from cffi import FFI
import _FDtomoC

class Mesh(object):
    def __init__(self, numberOfNode = None, igrid = None):
        self.numberOfNode = numberOfNode
        self.igrid = igrid

    @abc.abstractmethod
    def getSize(self):
        return NotImplemented

    @abc.abstractmethod
    def generateFineMesh(self):
        return NotImplemented

    @abc.abstractmethod
    def translate(self):
        return NotImplemented


class Mesh1D(Mesh):
    def getSize(self):
        return self.numberOfNode

    def translate(self):
        mesh = _FDtomoC.ffi.new("Mesh1D *", {'numberOfNode' : self.numberOfNode})
        igrid = _FDtomoC.ffi.new("int[]", self.igrid)
        mesh.igrid = igrid
        return mesh


class Mesh3D(Mesh):
    def getSize(self):
        size = 0
        for i in self.numberOfNode:
            size = size + i
        return size

    def translate(self):
        gridx = _FDtomoC.ffi.new("int[]", self.igrid[0])
        gridy = _FDtomoC.ffi.new("int[]", self.igrid[1])
        gridz = _FDtomoC.ffi.new("int[]", self.igrid[2])
        mesh = _FDtomoC.ffi.new("Mesh3D *", {'numberOfNode' : self.numberOfNode, 'gridx' : gridx, \
                                                'gridy' : gridy, 'gridz' : gridz})

        return mesh
