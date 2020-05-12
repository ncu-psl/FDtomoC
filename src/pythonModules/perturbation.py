import abc
from cffi import FFI
import _FDtomoC

class Perturbation(object):
    def __init__(self, perturbation = None):
        self.perturbationField = perturbation
        
