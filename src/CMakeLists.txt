SET(PROGS
	make1d
	c2f
	sphfd
	sphfdloc
	sprayderv
	runlsqr
	makenewmod
)

ADD_EXECUTABLE(FDtomo FDtomo.c)
FOREACH(prog ${PROGS})
	TARGET_LINK_LIBRARIES(FDtomo
		${prog}
	)
ENDFOREACH()

TARGET_INCLUDE_DIRECTORIES(FDtomo PUBLIC
	${MAIN_INCLUDE_DIR}
	${MAIN_INCLUDE_DIR}/FDtomo
)

