SET(PROGS
	make1d
	c2f
	sphfd
	sphfdloc
	sphrayderv
	runlsqr
	makenewmod
)

FOREACH(prog ${PROGS})
	ADD_EXECUTABLE(${prog}
		${prog}.c
	)
	TARGET_LINK_LIBRARIES(${prog}
		common_libraries
		m
	)
	TARGET_INCLUDE_DIRECTORIES(${prog} PUBLIC
		${MAIN_INCLUDE_DIR}
		${MAIN_INCLUDE_DIR}/common
	)
ENDFOREACH()

TARGET_LINK_LIBRARIES(sphrayderv
	sphrayderv_libraries
)

TARGET_LINK_LIBRARIES(runlsqr
	runlsqr_libraries
)
