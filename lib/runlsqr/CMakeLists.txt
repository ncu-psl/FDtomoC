ADD_LIBRARY(_runlsqr
	aprod.c
	makea.c
	normlz.c
	scopy.c
	snrm2.c
	sscal.c
)

include_directories(${MAIN_INCLUDE_DIR})
