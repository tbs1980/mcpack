# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 2.8)
project(TestHamiltonianMonteCarlo CXX)

find_package(Eigen3 3.1.2 REQUIRED)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

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

if(UNIT_TESTS)
    enable_testing()
endif()

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR})

add_executable(GaussLogPost GaussLogPost.cpp)
target_link_libraries(GaussLogPost ${Boost_LIBRARIES})
add_test(TestGaussLogPost GaussLogPost)

add_executable(KineticEnergy KineticEnergy.cpp)
target_link_libraries(KineticEnergy ${Boost_LIBRARIES})
add_test(TestKineticEnergy KineticEnergy)

add_executable(Integration Integration.cpp)
target_link_libraries(Integration ${Boost_LIBRARIES})
add_test(TestIntegration Integration)

add_executable(HMCProposal HMCProposal.cpp)
target_link_libraries(HMCProposal ${Boost_LIBRARIES})
add_test(TestHMCProposal HMCProposal)

add_executable(ClassicHMC ClassicHMC.cpp)
target_link_libraries(ClassicHMC ${Boost_LIBRARIES})
add_test(TestClassicHMC ClassicHMC)

add_executable(RuntimeControl RuntimeControl.cpp)
target_link_libraries(RuntimeControl ${Boost_LIBRARIES})
add_test(TestRuntimeControl RuntimeControl)

add_executable(IO IO.cpp)
target_link_libraries(IO ${Boost_LIBRARIES})
add_test(TestIO IO)

add_executable(Sampler Sampler.cpp)
target_link_libraries(Sampler ${Boost_LIBRARIES})
add_test(TestSampler Sampler)

add_subdirectory(MPI)

