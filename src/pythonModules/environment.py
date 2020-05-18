import abc
from cffi import FFI
import _FDtomoC

class CommonEnv(object):
    def __init__(self):
        self.common_env = None

    @abc.abstractmethod
    def create(self):
        return NotImplemented


class LocEnv(CommonEnv):
    def __init__(self, locEnvField = None, commonEnvField = None):
        self.locEnvField = locEnvField
        self.commonEnvField = commonEnvField

    def create(self, file = None):
        if (file != None):
            tmp = _FDtomoC.ffi.new("char[]", file.encode('ascii'))
            locEnvField = _FDtomoC.lib.setLocEnv(tmp)
            commonEnvField = _FDtomoC.lib.setCommonEnv(tmp)
            loc_env = LocEnv(locEnvField, commonEnvField)
            return loc_env


class SphraydervEnv(CommonEnv):
    def __init__(self, sphraydervEnvField = None, commonEnvField = None):
        self.sphraydervEnvField = sphraydervEnvField
        self.commonEnvField = commonEnvField

    def create(self, file = None):
        if (file != None):
            tmp = _FDtomoC.ffi.new("char[]", file.encode('ascii'))
            sphraydervEnvField = _FDtomoC.lib.setSphraydervEnv(tmp)
            commonEnvField = _FDtomoC.lib.setCommonEnv(tmp)
            sphrayderv_env = SphraydervEnv(sphraydervEnvField, commonEnvField)
            return sphrayderv_env


class RunlsqrEnv(CommonEnv):
    def __init__(self, runlsqrEnvField = None, commonEnvField = None):
        self.runlsqrEnvField = runlsqrEnvField
        self.commonEnvField = commonEnvField

    def create(self, file = None):
        if (file != None):
            runlsqr_env = RunlsqrEnv()
            tmp = _FDtomoC.ffi.new("char[]", file.encode('ascii'))
            runlsqrEnvField = _FDtomoC.lib.setRunlsqrEnv(tmp)
            commonEnvField = _FDtomoC.lib.setCommonEnv(tmp)
            runlsqr_env = RunlsqrEnv(runlsqrEnvField, commonEnvField)
            return runlsqr_env

class MakenewmodEnv(CommonEnv):
    def __init__(self, makenewmodEnvField = None, commonEnvField = None):
        self.makenewmodEnvField = makenewmodEnvField
        self.commonEnvField = commonEnvField

    def create(self, file = None):
        if (file != None):
            tmp = _FDtomoC.ffi.new("char[]", file.encode('ascii'))
            makenewmodEnvField = _FDtomoC.lib.setMakeNewmodEnv(tmp)
            commonEnvField = _FDtomoC.lib.setCommonEnv(tmp)
            makenewmod_env = MakenewmodEnv(makenewmodEnvField, commonEnvField)
            return makenewmod_env

