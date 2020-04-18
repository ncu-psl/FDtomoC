ADD_LIBRARY(_common
	parseprogs.c
	geographic_method.c
	string_process.c
	time_process.c
	# c_link_list.c
	earthquake_file_delimiter.c
	read_spec.c
	gridspec.c
	interpolation.c
	velocity_model.c
	grid.c
	station.c
)

include_directories(${MAIN_INCLUDE_DIR})
