#!/usr/bin/env python
from cffi import FFI
from subprocess import check_output
import os
from exportHeader import spec_header, grid_header, model_header, env_header, station_header, travelTimeTable
path = os.getcwd()
ffi=FFI()
header = spec_header.spec_header() + grid_header.grid_header() + model_header.model_header() \
        + env_header.env_header() + station_header.station_header() + travelTimeTable.travel_time_header()
ffi.cdef(header)
ffi.set_source("_FDtomoC",
     '''#include "FDtomo/make1d.h" 
        #include "FDtomo/sphfd.h"
        #include "FDtomo/runlsqr.h"
        #include "FDtomo/sphrayderv.h"
        #include "common/earthquake.h"
        #include "common/grid.h" 
        #include "common/vhead.h"
        #include "common/velocity_model.h"
        #include "common/station.h" 
    ''' ,
    #sources = [],
    include_dirs = [path + '/../../include'],
    libraries =['_common', 'runlsqr', 'sphrayderv', 'sphfdloc', 'sphfd'],
    library_dirs = [path + '/../../build/lib/common', path + '/../../build/lib/FDtomo']
    )



#if __name__ == "__main__":
ffi.compile(verbose=True)