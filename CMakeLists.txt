cmake_minimum_required(VERSION 3.14)

project(puma VERSION 0.1 LANGUAGES C CXX)

enable_language(Fortran)

find_package(VTKm REQUIRED)
find_package(ADIOS2 REQUIRED)
find_package(MPI REQUIRED)

add_subdirectory(src)
add_subdirectory(pumacli)

