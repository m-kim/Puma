add_executable(pumacli
	src/main.cxx
	)
set_target_properties(pumacli PROPERTIES CXX_STANDARD 14)
target_link_libraries(pumacli
        vtkm_cont
	puma
	)
