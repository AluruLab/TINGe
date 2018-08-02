cmake_minimum_required(VERSION 2.8.8)

set(PROJECT_NAME_STR tinge)
PROJECT(${PROJECT_NAME_STR} C CXX)
#TODO CHANGE THIS OPTION
set(CMAKE_VERBOSE_MAKEFILE ON)
# set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)

find_package(Threads REQUIRED)

# flags
set(CMAKE_CXX_FLAGS                "-Wall -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -std=c++11 -g")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -std=c++11")

# include jaz and mpix
include_directories("${PROJECT_SOURCE_DIR}/")
include_directories("${PROJECT_SOURCE_DIR}/src/")

#
# add MPI includes and libraries
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
message(STATUS "Found MPI:")
message(STATUS "    headers: ${MPI_INCLUDE_PATH}")
message(STATUS "    libs:    ${MPI_LIBRARIES}")
message(STATUS "    link flags: ${MPI_LINK_FLAGS}")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_LINK_FLAGS "${CMAKE_LINK_FLAGS} ${MPI_LINK_FLAGS}")

# Compile
set(PROJECT_EXEC ${PROJECT_NAME_STR}-mi)
ADD_EXECUTABLE(${PROJECT_EXEC} 
               src/main.cxx
               src/adj_build.cpp
               src/iomanip.cpp
               src/read_adj.cpp
               src/read_exp.cpp
               src/read_tfs.cpp
               src/transform_data.cpp
               src/write_adj.cpp
               src/write_lab.cpp)
# Link MPI library
target_link_libraries(${PROJECT_EXEC} "${MPI_LIBRARIES}")

# Compile
ADD_EXECUTABLE(residuals
               tools/residuals.cxx 
               src/iomanip.cpp
               src/read_exp.cpp)
target_link_libraries(residuals "${MPI_LIBRARIES}")


