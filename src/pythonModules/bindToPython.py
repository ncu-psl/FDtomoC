#!/usr/bin/env python
from cffi import FFI
from subprocess import check_output
import os
from pythonModules import spec_header, grid_header, model_header

ffi=FFI()
header = spec_header.spec_header() + grid_header.grid_header() + model_header.model_header()
ffi.cdef(header)
ffi.set_source("_test",
     '''#include "FDtomo/make1d.h" 
        #include "common/grid.h" 
        #include "common/vec/vec.h"
        #include "common/vhead.h"
        #include "common/velocity_model.h"
    ''' ,
    #sources = [],
    include_dirs = ['../include'],
    libraries =['_common'],
    library_dirs = ['../build/lib/common']
    )



#if __name__ == "__main__":
ffi.compile(verbose=True)