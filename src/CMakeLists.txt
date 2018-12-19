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
	parseprogs.c
	string_process.c
	g_method.c

	# sphfdloc sphrayderv
	time_process.c
	earthquake_file_delimiter.c

	#runlsqr
	makea.c
	scopy.c
	aprod.c
	normlz.c
	sscal.c
	snrm2.c
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
