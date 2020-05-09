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
    def __init__(self):
        self.loc_env = None

    def create(self, file = None):
        if (file != None):
            loc_env = LocEnv()
            tmp = _FDtomoC.ffi.new("char[]", file)
            loc_env.loc_env = _FDtomoC.lib.setLocEnv(tmp)
            loc_env.common_env = _FDtomoC.lib.setCommonEnv(tmp)
            return loc_env


class SphraydervEnv(CommonEnv):
    def __init__(self):
        self.sphrayderv_env = None

    def create(self, file = None):
        if (file != None):
            sphrayderv_env = SphraydervEnv()
            tmp = _FDtomoC.ffi.new("char[]", file)
            sphrayderv_env.sphrayderv_env = _FDtomoC.lib.setSphraydervEnv(tmp)
            sphrayderv_env.common_env = _FDtomoC.lib.setCommonEnv(tmp)
            return sphrayderv_env


class RunlsqrEnv(CommonEnv):
    def __init__(self):
        self.runlsqr_env = None

    def create(self, file = None):
        if (file != None):
            runlsqr_env = RunlsqrEnv()
            tmp = _FDtomoC.ffi.new("char[]", file)
            runlsqr_env.runlsqr_env = _FDtomoC.lib.setRunlsqrEnv(tmp)
            runlsqr_env.common_env = _FDtomoC.lib.setCommonEnv(tmp)
            return runlsqr_env

