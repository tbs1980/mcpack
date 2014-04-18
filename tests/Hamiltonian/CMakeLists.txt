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

find_package(Eigen3 3.1.2 REQUIRED)
find_package(Boost REQUIRED COMPONENTS unit_test_framework)

if(UNIX)
	if(CMAKE_COMPILER_IS_GNUCXX)
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -pedantic -Wall -Wextra -Wfatal-errors -g")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -O3 -march=native")
	endif(CMAKE_COMPILER_IS_GNUCXX)
endif(UNIX)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR})

add_executable(KineticEnergy KineticEnergy.cpp)
target_link_libraries(KineticEnergy ${Boost_LIBRARIES})
add_test(TestKineticEnergy GaussKineticEnergy)

