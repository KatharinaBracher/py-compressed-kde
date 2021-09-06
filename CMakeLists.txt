cmake_minimum_required(VERSION 2.8)

ENABLE_LANGUAGE(CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wunused -std=c++17")

option(BUILD_SHARED_LIBS "Build shared library" OFF)

project(DECODER)

include(FetchContent)
FetchContent_Declare(
		yaml-cpp
		GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
		GIT_SHALLOW	ON
		GIT_TAG yaml-cpp-0.6.3
    )
    
set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "disable yaml tests")
message(STATUS "Populating yaml-cpp...")
FetchContent_MakeAvailable(yaml-cpp)
message(STATUS "Done.")
find_package(HDF5 REQUIRED)

include_directories( ${HDF5_INCLUDE_DIRS} )
include_directories( ${YAML_CPP_INCLUDE_DIR} )
include_directories( src )
include_directories( ext/HighFive-1.4/include )

file(GLOB sources "src/*.cpp")
file(GLOB header_files "src/*.hpp")

add_library(compressed_decoder ${sources})
target_link_libraries( compressed_decoder yaml-cpp ${HDF5_C_LIBRARIES})

set(INCLUDE_INSTALL_ROOT_DIR ${CMAKE_INSTALL_PREFIX}/include)
set(INCLUDE_INSTALL_DIR ${INCLUDE_INSTALL_ROOT_DIR}/compressed_decoder)

install (TARGETS compressed_decoder DESTINATION lib)
install (FILES ${header_files} DESTINATION ${INCLUDE_INSTALL_DIR})
