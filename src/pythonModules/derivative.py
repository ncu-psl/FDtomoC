import abc
from cffi import FFI
import _FDtomoC

class Derivative(object):
    def __init__(self, derivative = None):
        self.derivativeField = derivative
        
