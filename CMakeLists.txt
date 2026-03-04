cmake_minimum_required(VERSION 3.16)

project(FEM)

set(CMAKE_CXX_STANDARD 20)

find_package(OpenMP REQUIRED)

include_directories(include)

add_library(fem_libs
    src/funciones_1d.cpp
    src/utils.cpp
)

# Crear ejecutable
add_executable(main_exec main/fem.cpp)

#juntar librería
target_link_libraries(main_exec fem_libs)

#juntar librería de openmp
target_link_libraries(main_exec
    fem_libs
    OpenMP::OpenMP_CXX
)