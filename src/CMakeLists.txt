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
		_common
	)
	IF(UNIX)
		TARGET_LINK_LIBRARIES(${prog}
			m
		)
	ENDIF(UNIX)
	TARGET_INCLUDE_DIRECTORIES(${prog} PUBLIC
		${MAIN_INCLUDE_DIR}
		${MAIN_INCLUDE_DIR}/common
	)
ENDFOREACH()

TARGET_LINK_LIBRARIES(sphrayderv
	_sphrayderv
)

TARGET_LINK_LIBRARIES(runlsqr
	_runlsqr
)
