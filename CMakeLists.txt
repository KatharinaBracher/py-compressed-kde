cmake_minimum_required(VERSION 2.8)

ENABLE_LANGUAGE(CXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wunused -std=c++17")

option(BUILD_SHARED_LIBS "Build shared library" OFF)

project(DECODER)

include(FetchContent)

FetchContent_Declare(
    flatbuffers
    GIT_REPOSITORY "https://github.com/google/flatbuffers"
    GIT_TAG        v2.0.0
)

FetchContent_Declare(
		yaml-cpp
		GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
		GIT_SHALLOW	ON
		GIT_TAG yaml-cpp-0.6.3
    )

set(FLATBUFFERS_BUILD_FLATC ON CACHE INTERNAL "Enable the build of the flatbuffers compiler")
set(FLATBUFFERS_BUILD_TESTS OFF CACHE INTERNAL "Enable tests from flatbuffers lib.")
set(FLATBUFFERS_BUILD_FLATHASH OFF CACHE INTERNAL "Enable the build of flathash")
set(FLATBUFFERS_ENABLE_PCH OFF CACHE INTERNAL
		"Enable precompile headers support for 'flatbuffers' and 'flatc'. \"
		Only work if CMake support 'target_precompile_headers'. \"
		Only work if CMake supports 'target_precompile_headers'. \"
		This can speed up compilation time."
		)

FetchContent_MakeAvailable(flatbuffers)

include_directories(BEFORE SYSTEM "${flatbuffers_SOURCE_DIR}/include")

set(YAML_CPP_BUILD_TESTS OFF CACHE BOOL "disable yaml tests")
message(STATUS "Populating yaml-cpp...")
FetchContent_MakeAvailable(yaml-cpp)
message(STATUS "Done.")
find_package(HDF5 REQUIRED)

include_directories( ${HDF5_INCLUDE_DIRS} )
include_directories( ${YAML_CPP_INCLUDE_DIR} )
include_directories( src )
include_directories( ext/HighFive-1.4/include )
include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/datatype_generated.h
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/datatype.fbs
    COMMAND flatc --cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/datatype.fbs
)

add_custom_target(datatypebuffer DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/datatype_generated.h)

file(GLOB sources "src/*.cpp")
file(GLOB header_files "src/*.hpp")

add_library(compressed_decoder ${sources})
add_dependencies(compressed_decoder datatypebuffer)
target_link_libraries( compressed_decoder yaml-cpp flatbuffers ${HDF5_C_LIBRARIES})

set(INCLUDE_INSTALL_ROOT_DIR ${CMAKE_INSTALL_PREFIX}/include)
set(INCLUDE_INSTALL_DIR ${INCLUDE_INSTALL_ROOT_DIR}/compressed_decoder)

install (TARGETS compressed_decoder DESTINATION lib)
install (FILES ${header_files} DESTINATION ${INCLUDE_INSTALL_DIR})
