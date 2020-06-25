set(HEADER_LIST ${puma_SOURCE_DIR}/include/XgcExtrudeCompute.h
${puma_SOURCE_DIR}/include/XgcExtrudeMesh.h
${puma_SOURCE_DIR}/include/XgcExtrude.h
${puma_SOURCE_DIR}/include/TurbulenceWorklets.h
)

set(SOURCE_LIST
	main.cxx
        )
add_executable(pumacli ${HEADER_LIST} ${SOURCE_LIST}
	)
target_include_directories(pumacli PUBLIC "../include")
set_target_properties(pumacli PROPERTIES CXX_STANDARD 14)
set(LIBRARY_LIST 
puma
vtkm_cont
    vtkm_rendering
    kittie::kittie
    MPI::MPI_C
    MPI::MPI_CXX
    MPI::MPI_Fortran
    adios2::adios2
	)
set_target_properties(pumacli PROPERTIES CXX_STANDARD 14)
if(VTKm_ENABLE_CUDA)
  set_source_files_properties(${SOURCE_LIST} PROPERTIES LANGUAGE "CUDA")
  set_target_properties(pumacli PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  set_target_properties(pumacli PROPERTIES
    CUDA_STANDARD 14)
  set(LIBRARY_LIST ${LIBRARY_LIST} vtkm::cuda)
endif(VTKm_ENABLE_CUDA)
target_link_libraries(pumacli
        ${LIBRARY_LIST})
install(
TARGETS pumacli
RUNTIME DESTINATION bin
)
