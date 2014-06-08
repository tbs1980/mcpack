# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 2.8)
project(TestHamiltonianMonteCarlo_MPI CXX)

find_package(MPI)
find_package(Boost COMPONENTS unit_test_framework mpi serialization)

if(UNIT_TESTS)
    enable_testing()
endif()


if(MPI_CXX_FOUND)
    if(Boost_FOUND)
        include_directories(${MPI_CXX_INCLUDE_PATH})
        include_directories(${Boost_INCLUDE_DIRS})
        include_directories(${CMAKE_SOURCE_DIR})

        add_executable(MPI_Sampler MPI_Sampler.cpp)
        target_link_libraries(MPI_Sampler ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})
        add_test(TestMPI_Sampler MPI_Sampler)
    endif()
endif()
