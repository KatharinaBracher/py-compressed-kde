cmake_minimum_required(VERSION 2.8)

ENABLE_LANGUAGE(CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wunused -std=c++11")

option(BUILD_SHARED_LIBS "Build shared library" ON)

project(DECODER)

find_package( YAML-CPP REQUIRED )

include_directories( ${YAML_CPP_INCLUDE_DIR} )
include_directories( 'src' )

file(GLOB sources "src/*.cpp")
file(GLOB header_files "src/*.hpp")

add_library(compressed_decoder ${sources})
target_link_libraries( compressed_decoder yaml-cpp )


set(INCLUDE_INSTALL_ROOT_DIR ${CMAKE_INSTALL_PREFIX}/include)
set(INCLUDE_INSTALL_DIR ${INCLUDE_INSTALL_ROOT_DIR}/compressed_decoder)

install (TARGETS compressed_decoder DESTINATION lib)
install (FILES ${header_files} DESTINATION ${INCLUDE_INSTALL_DIR})
