ADD_LIBRARY(_runlsqr STATIC
	aprod.c
	makea.c
	normlz.c
	scopy.c
	snrm2.c
	sscal.c
)

include_directories(${MAIN_INCLUDE_DIR})
