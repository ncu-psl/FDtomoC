cmake_minimum_required(VERSION 2.8)
project(FDtomo)
include_directories("../include")

SET(CMAKE_CXX_STANDRAD 14)

if(CMAKE_C_COMPILER_ID MATCHES "GNU")
	SET(CMAKE_C_FLAGS " -O3 -Wall -std=c99 -fopenmp -mcmodel=large -m64 -lm -fdata-sections -ffunction-sections -Wl,--gc-sections")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
	SET(CMAKE_C_FLAGS " -O3 -Wall -std=c99 -qopenmp -mcmodel=large -m64 -ipo -align")
endif()

SET(PROGS
	make1d c2f sphfdloc sphrayderv runlsqr makenewmod
)

SET(PROG_LIB
	# common lib
	c_parseprogs.c
	c_string_process.c
	c_g_method.c

	# sphfdloc sphrayderv
	c_time_process.c
	c_earthquake_file_delimiter.c

	#runlsqr
	c_makea.c
	c_scopy.c
	c_aprod.c
	c_normlz.c
	c_sscal.c
	c_snrm2.c
)

foreach(prog ${PROGS})
	add_executable(${prog}
		${prog}.c
		${PROG_LIB}
	)
	target_link_libraries(${prog} m)
endforeach()

add_executable(sphfd
	sphfd.c
)
