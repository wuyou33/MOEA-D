cmake_minimum_required(VERSION 3.8)
project(MOEA_D)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp loader.hpp DataSet DataSet/portreb1.txt moead.h)
add_executable(MOEA_D ${SOURCE_FILES} loader.hpp)