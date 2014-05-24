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
