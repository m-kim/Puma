cmake_minimum_required(VERSION 3.14)

project(TestHarness VERSION 0.1 LANGUAGES C CXX)

enable_language(Fortran)

find_package(VTKm REQUIRED)
find_package(ADIOS2 REQUIRED)
find_package(MPI REQUIRED)
set(TestHarnessSrc
        main.cxx)
add_executable(TestHarness ${TestHarnessSrc})
set_target_properties(TestHarness PROPERTIES CXX_STANDARD 14)

target_link_libraries(TestHarness PUBLIC
    vtkm_rendering
    MPI::MPI_C
    MPI::MPI_CXX
    MPI::MPI_Fortran
    adios2::adios2 )
