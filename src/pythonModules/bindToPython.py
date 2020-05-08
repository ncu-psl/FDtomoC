#!/usr/bin/env python
from cffi import FFI
from subprocess import check_output
import os
from exportHeader import spec_header, grid_header, model_header
path = os.getcwd()
ffi=FFI()
header = spec_header.spec_header() + grid_header.grid_header() + model_header.model_header()
ffi.cdef(header)
ffi.set_source("_FDtomoC",
     '''#include "FDtomo/make1d.h" 
        #include "common/grid.h" 
        #include "common/vec/vec.h"
        #include "common/vhead.h"
        #include "common/velocity_model.h"
    ''' ,
    #sources = [],
    include_dirs = [path + '/../../include'],
    libraries =['_common'],
    library_dirs = [path + '/../../build/lib/common']
    )



#if __name__ == "__main__":
ffi.compile(verbose=True)