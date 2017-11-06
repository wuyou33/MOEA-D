cmake_minimum_required(VERSION 3.5)
project(MOEA_D)

set(CMAKE_CXX_STANDARD 17)
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRES} ${GSLCBLAS_INCLUDE_DIRS})
link_libraries(${GSL_LIBRARIES} ${GSLBLAS_LIBRARIES})



set(SOURCE_FILES main.cpp loader.hpp DataSet DataSet/portreb1.txt moead.h randG.h)
add_executable(MOEA_D ${SOURCE_FILES} loader.hpp randG.h)
target_link_libraries (MOEA_D ${GSL_LIBRARIES} ${GSLBLAS_LIBRARIES})