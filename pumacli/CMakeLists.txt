add_executable(pumacli
	main.cxx
	)
set_target_properties(pumacli PROPERTIES CXX_STANDARD 14)
target_link_libraries(pumacli
	puma
	)
