# 
# 
# Copyright (C) 2014 Sreekumar Thaithara Balan
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#


cmake_minimum_required(VERSION 2.8)
project(TestHamiltonianMonteCarlo CXX)

find_package(MPI REQUIRED)
find_package(Eigen3 3.1.2 REQUIRED)
find_package(Boost REQUIRED COMPONENTS unit_test_framework mpi serialization)

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

include_directories(${MPI_CXX_INCLUDE_PATH})
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

add_executable(MPI_Sampler MPI_Sampler.cpp)
target_link_libraries(MPI_Sampler ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})
add_test(TestMPI_Sampler MPI_Sampler)
