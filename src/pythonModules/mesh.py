import abc
from cffi import FFI
import _FDtomoC

class Mesh(object):
    def __init__(self, numberOfNode = None, igrid = None):
        self.numberOfNode = numberOfNode
        self.igrid = igrid
        self.meshField = None

    @abc.abstractmethod
    def create(self):
        return NotImplemented

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
    def create(self, file = None, numberOfNode = None, igrid = None):
        mesh = Mesh1D(numberOfNode, igrid)
        tmp = _FDtomoC.ffi.new("int[]", igrid)
        mesh.meshField = _FDtomoC.lib.createMesh1D(numberOfNode, tmp)
        return mesh

    def getSize(self):
        return self.numberOfNode

    def translate(self):
        mesh = _FDtomoC.ffi.new("Mesh1D *", {'numberOfNode' : self.numberOfNode})
        igrid = _FDtomoC.ffi.new("int[]", self.igrid)
        mesh.igrid = igrid
        self.mesh = mesh


class Mesh3D(Mesh):
    def create(self, file = None, numberOfNode = None, igrid = None):
        if (file != None):
            mesh = Mesh3D()
            tmp = _FDtomoC.ffi.new("char[]", file)
            mesh.meshField = _FDtomoC.lib.setMesh3D(tmp)
            return mesh


    def getSize(self):
        size = 0
        for i in self.numberOfNode:
            size = size + i
        return size

    def generateFineMesh(self):
        mesh = Mesh3D()
        mesh.meshField = _FDtomoC.lib.generateFineMesh(self.meshField)
        return mesh

    def translate(self):
        gridx = _FDtomoC.ffi.new("int[]", self.igrid[0])
        gridy = _FDtomoC.ffi.new("int[]", self.igrid[1])
        gridz = _FDtomoC.ffi.new("int[]", self.igrid[2])
        mesh = _FDtomoC.ffi.new("Mesh3D *", {'numberOfNode' : self.numberOfNode, 'gridx' : gridx, \
                                                'gridy' : gridy, 'gridz' : gridz})

        self.mesh = mesh 
