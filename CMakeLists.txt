cmake_minimum_required(VERSION 3.14)

project(Puma VERSION 0.1 LANGUAGES C CXX)

enable_language(Fortran)

find_package(VTKm REQUIRED)
find_package(ADIOS2 REQUIRED)
find_package(MPI REQUIRED)
set(PumaSrc
        XgcExtrudeCompute.cxx
        XgcExtrudeMesh.cxx
        main.cxx)
add_executable(Puma ${PumaSrc})
set_target_properties(Puma PROPERTIES CXX_STANDARD 14)

target_link_libraries(Puma PUBLIC
    vtkm_cont
    vtkm_rendering
    MPI::MPI_C
    MPI::MPI_CXX
    MPI::MPI_Fortran
    adios2::adios2 )
