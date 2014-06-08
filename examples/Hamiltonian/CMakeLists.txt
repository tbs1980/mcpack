# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


cmake_minimum_required(VERSION 2.8)
project(HMCExamples CXX)

find_package(Eigen3 3.1.2 REQUIRED)
find_package(Boost REQUIRED)

if(UNIX)
    if(CMAKE_COMPILER_IS_GNUCXX)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -pedantic -Wall -Wextra -Wfatal-errors -g")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -O3 -march=native")
    endif(CMAKE_COMPILER_IS_GNUCXX)
    
    #is the compiler clang?
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -pedantic -Wall -Wextra -g")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -O3 -march=native")
    elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=gnu++11 -pedantic -w3 -g -I/usr/include/x86_64-linux-gnu/c++/4.8")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -I/usr/include/x86_64-linux-gnu/c++/4.8 -std=gnu++98 -O3 -march=native ")
    endif()     
endif(UNIX)

include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR})

add_executable(Rosenbrock Rosenbrock.cpp)

