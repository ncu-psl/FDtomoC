SET(PROGS
	make1d
	c2f
	sphfd
	sphfdloc
	sphrayderv
	runlsqr
	makenewmod
)

ADD_EXECUTABLE(Fdtomo FDtomo.c)
FOREACH(prog ${PROGS})
	TARGET_LINK_LIBRARIES(Fdtomo
		${prog}
	)
ENDFOREACH()

TARGET_INCLUDE_DIRECTORIES(Fdtomo PUBLIC
	${MAIN_INCLUDE_DIR}
	${MAIN_INCLUDE_DIR}/FDtomo
)
