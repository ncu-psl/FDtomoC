import abc
from cffi import FFI
import _FDtomoC

class VelocityModel(object):
    def __init__(self, coordinate = None, velocity = None):
        self.coordinate = coordinate
        self.velocity = velocity
        self.modelField = None
        self.modelFieldPtr = None

    @abc.abstractmethod
    def create(self):
        return NotImplemented

    @abc.abstractmethod
    def transform(self):
        return NotImplemented


class VelocityModel1D(VelocityModel):
    def setVelocityModel(self, file, vpModel, vsModel):
        interp = _FDtomoC.ffi.new("char *")
        vpModelField = _FDtomoC.ffi.new("velocityModel1D *")
        vsModelField = _FDtomoC.ffi.new("velocityModel1D *")
        tmp = _FDtomoC.ffi.new("char[]", file.encode('ascii'))
        _FDtomoC.lib.readVelocityModel1D(tmp, vpModelField, vsModelField, interp)
        vpModel.modelField = vpModelField[0]
        vsModel.modelField = vsModelField[0]

    def transform(self, coordinate1D):
        model = VelocityModel1D()
        tmp = ""
        for i in range(coordinate1D.coordinateField.mesh.numberOfNode):
            tmp = tmp + 'I'
        tmp = _FDtomoC.ffi.new("char[]", tmp.encode('ascii'))
        model.modelField = _FDtomoC.lib.transform1D(coordinate1D.coordinateField, self.modelField, tmp)
        return model        


class VelocityModel3D(VelocityModel):
    def create(self, coordinate3D, model1D):
        model3D = VelocityModel3D()
        model3D.modelField = _FDtomoC.lib.create3DModel(coordinate3D.coordinateField, model1D.modelField)
        return model3D

    def transform(self, coordinate3D):
        model3D = VelocityModel3D()
        model3D.modelField = _FDtomoC.lib.transform3D(coordinate3D.coordinateField, self.modelField)
        return model3D

    def makeNewModel(self, coordinate, vp_model, vs_model, perturbation, table_size, new_model_env):
        corField = coordinate.coordinateField

        vpModelFieldPtr = _FDtomoC.ffi.new("velocityModel3D *", vp_model.modelField)
        vsModelFieldPtr = _FDtomoC.ffi.new("velocityModel3D *", vs_model.modelField)

        perturbationField = perturbation.perturbationField

        makenewmodEnvField = new_model_env.makenewmodEnvField
        commonEnvField = new_model_env.commonEnvField

        _FDtomoC.lib.makenewmod(corField, vpModelFieldPtr, vsModelFieldPtr, perturbationField, table_size, makenewmodEnvField, commonEnvField)
        vp_model.modelField = vpModelFieldPtr[0]
        vs_model.modelField = vsModelFieldPtr[0]
