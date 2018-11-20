cmake_minimum_required(VERSION 2.8)
project(FDtomo)
include_directories("../include")

SET(CMAKE_CXX_STANDRAD 14)

if(CMAKE_C_COMPILER_ID MATCHES "GNU")
	SET(CMAKE_C_FLAGS " -O3 -Wall -std=c99 -fopenmp -mcmodel=large -m64 -lm -fdata-sections -ffunction-sections -Wl,--gc-sections")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
	SET(CMAKE_C_FLAGS " -O3 -Wall -std=c99 -qopenmp -mcmodel=large -m64 -ipo -align")
endif()

add_executable(make1d
	make1d.c
	c_parseprogs.c
	c_g_method.c
	c_string_process.c
)
target_link_libraries(make1d m)

add_executable(c2f
	c2f.c
	c_parseprogs.c
	c_string_process.c
)
target_link_libraries(c2f m)

add_executable(sphfdloc
	sphfdloc.c
	c_parseprogs.c
	c_g_method.c
	c_time_process.c
	c_earthquake_file_delimiter.c
	c_string_process.c
)
target_link_libraries(sphfdloc m)

add_executable(runlsqr
	runlsqr.c
	c_parseprogs.c
	c_makea.c
	c_scopy.c
	c_aprod.c
	c_normlz.c
	c_sscal.c
	c_snrm2.c
)
target_link_libraries(runlsqr m)

add_executable(makenewmod
	makenewmod.c
	c_parseprogs.c
	c_string_process.c
)
target_link_libraries(makenewmod m)
