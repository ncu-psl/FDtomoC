import abc
from cffi import FFI
from point import Point
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
    def getField(self):
        return NotImplemented

    @abc.abstractmethod
    def getClass(self):
        return NotImplemented

class Mesh1D(Mesh):
    def create(self, file = None, numberOfNode = None, igrid = None):
        mesh = Mesh1D(numberOfNode, igrid)
        mesh.meshField = mesh.getField()
        return mesh

    def getSize(self):
        return self.numberOfNode

    def getField(self):
        igrid = _FDtomoC.ffi.new("int[]", self.igrid)
        #meshFieldPtr = _FDtomoC.ffi.new("Mesh1D *", {'numberOfNode' : self.numberOfNode, 'igrid' : igrid})
        meshField = _FDtomoC.lib.createMesh1D(self.numberOfNode, igrid)
        return meshField

    def getClass(self):
        self.numberOfNode = self.meshField.numberOfNode
        self.igrid = _FDtomoC.ffi.unpack(self.meshField.igrid)

class Mesh3D(Mesh):
    def __init__(self, numberOfNode = None, gridx = None, gridy = None, gridz = None):
        if numberOfNode == None:
            self.numberOfNode = Point()
        else:
            self.numberOfNode = numberOfNode
        self.gridx = gridx
        self.gridy = gridy
        self.gridz = gridz

    def create(self, file = None, numberOfNode = None, gridx = None, gridy = None, gridz = None):
        if (file != None):
            mesh = Mesh3D()
            tmp = _FDtomoC.ffi.new("char[]", file.encode('ascii'))
            mesh.meshField = _FDtomoC.lib.setMesh3D(tmp)
            mesh.numberOfNode.pointField = mesh.meshField.numberOfNode
            mesh.getClass()
            return mesh

        mesh = Mesh3D(numberOfNode, gridx, gridy, gridz)
        mesh.meshField = mesh.getField()
        mesh.numberOfNode.pointField = mesh.meshField.numberOfNode
        return mesh

    def getSize(self):
        size = 0
        for i in self.numberOfNode:
            size = size + i
        return size

    def generateFineMesh(self):
        mesh = Mesh3D()
        mesh.meshField = _FDtomoC.lib.generateFineMesh(self.meshField)
        mesh.numberOfNode.pointField = mesh.meshField.numberOfNode
        mesh.getClass()
        return mesh

    def getField(self):
        pointField = self.numberOfNode.getField()
        gridx = _FDtomoC.ffi.new("int[]", self.gridx)
        gridy = _FDtomoC.ffi.new("int[]", self.gridy)
        gridz = _FDtomoC.ffi.new("int[]", self.gridz)
        meshField = _FDtomoC.lib.createMesh3D(pointField, gridx, gridy, gridz)
        '''
        meshFieldPtr = _FDtomoC.ffi.new("Mesh3D *", {'numberOfNode' : pointField, 'gridx' : gridx, \
                                                'gridy' : gridy, 'gridz' : gridz})
        '''
        return meshField

    def getClass(self):
        self.gridx = _FDtomoC.ffi.unpack(self.meshField.gridx, int(self.meshField.numberOfNode.x) - 1)
        self.gridy = _FDtomoC.ffi.unpack(self.meshField.gridy, int(self.meshField.numberOfNode.y) - 1)
        self.gridz = _FDtomoC.ffi.unpack(self.meshField.gridz, int(self.meshField.numberOfNode.z) - 1)
        self.numberOfNode.getClass()
