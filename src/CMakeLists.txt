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
target_include_directories(puma PUBLIC ../include)
set_target_properties(puma PROPERTIES CXX_STANDARD 14)

target_link_libraries(puma PUBLIC
    vtkm_cont
    vtkm_rendering
    MPI::MPI_C
    MPI::MPI_CXX
    MPI::MPI_Fortran
    adios2::adios2 )
