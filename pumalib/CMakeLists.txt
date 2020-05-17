cmake_minimum_required(VERSION 3.14)

project(PumaLib VERSION 0.1 LANGUAGES C CXX)


set(PumaSrc
        src/XgcExtrudeCompute.cxx
        src/XgcExtrudeMesh.cxx
	)
add_library(puma ${PumaSrc})

target_include_directories(puma PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
	$<INSTALL_INTERFACE:include>
	PRIVATE src
	)
set_target_properties(puma PROPERTIES CXX_STANDARD 14)

target_link_libraries(puma PUBLIC
    vtkm_cont
    vtkm_rendering
    MPI::MPI_C
    MPI::MPI_CXX
    MPI::MPI_Fortran
    adios2::adios2 )


