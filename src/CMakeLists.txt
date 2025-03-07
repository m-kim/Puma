set(HEADER_LIST ${puma_SOURCE_DIR}/include/XgcExtrudeCompute.h
${puma_SOURCE_DIR}/include/XgcExtrudeMesh.h
${puma_SOURCE_DIR}/include/XgcExtrude.h
${puma_SOURCE_DIR}/include/TurbulenceWorklets.h
)

set(SOURCE_LIST
        XgcExtrude.cxx
        XgcExtrudeCompute.cxx
        XgcExtrudeMesh.cxx
        )

add_library(puma ${SOURCE_LIST} ${HEADER_LIST})
target_include_directories(puma PUBLIC ../include )
set_target_properties(puma PROPERTIES CXX_STANDARD 14)

target_link_libraries(puma PUBLIC
    vtkm_cont
    vtkm_rendering
    kittie::kittie
    MPI::MPI_C
    MPI::MPI_CXX
    MPI::MPI_Fortran
    adios2::adios2 )

if (VTKm_ENABLE_CUDA)
  set_source_files_properties(${SOURCE_LIST} PROPERTIES LANGUAGE "CUDA")
  set_target_properties(puma PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  set_target_properties(puma PROPERTIES CUDA_STANDARD 14)

  target_link_libraries(puma PUBLIC vtkm::cuda)
    
endif(VTKm_ENABLE_CUDA)

