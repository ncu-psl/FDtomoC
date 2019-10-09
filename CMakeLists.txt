CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(Fdtomo)

SET(CMAKE_CXX_STANDRAD 14)
SET(MAIN_SRC_DIR     ${CMAKE_CURRENT_SOURCE_DIR})
SET(MAIN_INCLUDE_DIR ${MAIN_SRC_DIR}/include)
SET(MAIN_LIB_DIR     ${MAIN_SRC_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# set default build type = Release
if(NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release)
endif()

# gcc flag
SET(gcc_flags " -Wall -std=c99 -fopenmp -mcmodel=large -m64 -lm")
SET(gcc_flags_debug "${gcc_flags} -g")
SET(gcc_flags_release "${gcc_flags} -O3 -fdata-sections -ffunction-sections -Wl,--gc-sections")

# icc flag
SET(icc_flags " -Wall -std=c99 -qopenmp -mcmodel=large -m64 -align")
SET(icc_flags_debug "${icc_flags} -g")
SET(icc_flags_release "${icc_flags} -O3")

# if using dynamic linking, enable IPO(InterProcedural Optimization)
if(BUILD_SHARED_LIBS)
	SET(INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# determine compiler
if(CMAKE_C_COMPILER_ID MATCHES "GNU")
	SET(CMAKE_C_FLAGS ${gcc_flags})
	SET(CMAKE_C_FLAGS_DEBUG ${gcc_flags_debug})
	SET(CMAKE_C_FLAGS_RELEASE ${gcc_flags_release})
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
	SET(CMAKE_C_FLAGS ${icc_flags})
	SET(CMAKE_C_FLAGS_DEBUG ${icc_flags_debug})
	SET(CMAKE_C_FLAGS_RELEASE ${icc_flags_release})
endif()

ADD_SUBDIRECTORY(${MAIN_INCLUDE_DIR})
ADD_SUBDIRECTORY(${MAIN_LIB_DIR})
ADD_SUBDIRECTORY(src)
