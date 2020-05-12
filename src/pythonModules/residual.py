import abc
from cffi import FFI
import _FDtomoC

class ResidualVector(object):
    def __init__(self, residual_vector = None):
        self.residualField = residual_vector

        