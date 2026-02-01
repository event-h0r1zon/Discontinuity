cmake_minimum_required(VERSION 3.20)

project(riemann_solvers VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(riemann
    src/main.cpp
    src/riemann.cpp
    src/solver.cpp
)

target_include_directories(riemann PRIVATE src)