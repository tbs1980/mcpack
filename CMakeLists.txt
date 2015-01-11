# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 2.8)
project(mcpack)

#add the customised package searches to the module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# C and C++ compiler flags
include(CheckCXXCompilerFlag)
include(CheckCCompilerFlag)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(NOT cmake_build_type_tolower STREQUAL "debug" AND NOT cmake_build_type_tolower STREQUAL "release")
    message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug and Release.")
endif()

find_package(StandardMathLibrary)

set(MCPACKL_TEST_CUSTOM_LINKER_FLAGS  "" CACHE STRING "Additional linker flags when linking unit tests.")
set(MCPACKL_TEST_CUSTOM_CXX_FLAGS     "" CACHE STRING "Additional compiler flags when compiling unit tests.")

set(MCPACKL_STANDARD_LIBRARIES_TO_LINK_TO "")

if(NOT STANDARD_MATH_LIBRARY_FOUND)
    message(FATAL_ERROR
    "Can't link to the standard math library. Please report to the Eigen developers, telling them about your platform.")

else()
    if(MCPACKL_STANDARD_LIBRARIES_TO_LINK_TO)
        set(MCPACKL_STANDARD_LIBRARIES_TO_LINK_TO "${MCPACKL_STANDARD_LIBRARIES_TO_LINK_TO} ${STANDARD_MATH_LIBRARY}")
    else()
        set(MCPACKL_STANDARD_LIBRARIES_TO_LINK_TO "${STANDARD_MATH_LIBRARY}")
    endif()
endif()

if(MCPACKL_STANDARD_LIBRARIES_TO_LINK_TO)
    message(STATUS "Standard libraries to link to explicitly: ${MCPACKL_STANDARD_LIBRARIES_TO_LINK_TO}")
else()
    message(STATUS "Standard libraries to link to explicitly: none")
endif()

#add the customised package searches to the module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

macro(ei_add_cxx_compiler_flag FLAG)
    string(REGEX REPLACE "-" "" SFLAG1 ${FLAG})
    string(REGEX REPLACE "\\+" "p" SFLAG ${SFLAG1})
    check_cxx_compiler_flag(${FLAG} COMPILER_SUPPORT_${SFLAG})
    if(COMPILER_SUPPORT_${SFLAG})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FLAG}")
    endif()
endmacro(ei_add_cxx_compiler_flag)

macro(ei_add_c_compiler_flag FLAG)
    string(REGEX REPLACE "-" "" SFLAG1 ${FLAG})
    string(REGEX REPLACE "\\+" "p" SFLAG ${SFLAG1})
    check_c_compiler_flag(${FLAG} COMPILER_SUPPORT_${SFLAG})
    if(COMPILER_SUPPORT_${SFLAG})
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${FLAG}")
    endif()
endmacro(ei_add_c_compiler_flag)


if(NOT MSVC)
    # We assume that other compilers are partly compatible with GNUCC
    set(CMAKE_CXX_FLAGS_DEBUG "-g3")
    set(CMAKE_CXX_FLAGS_RELEASE "-g0 -O3")

    set(CMAKE_C_FLAGS_DEBUG "-g3")
    set(CMAKE_C_FLAGS_RELEASE "-g0 -O3")

    ei_add_cxx_compiler_flag("-pedantic")
    ei_add_cxx_compiler_flag("-Wall")
    ei_add_cxx_compiler_flag("-Wextra")
    ei_add_cxx_compiler_flag("-std=c++11")

    ei_add_c_compiler_flag("-pedantic")
    ei_add_c_compiler_flag("-Wall")
    ei_add_c_compiler_flag("-Wextra")
    ei_add_c_compiler_flag("-std=c11")
endif(NOT MSVC)


find_package(Eigen3 3.1.2 REQUIRED)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR})


enable_testing()
add_subdirectory(tests)
