# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 2.8)
project(TestHamiltonianMonteCarlo CXX)

add_subdirectory(utils)

add_executable(Test_kineticEnergy kineticEnergy.cpp)
target_link_libraries(Test_kineticEnergy ${Boost_LIBRARIES})
add_test(TestKineticEnergy Test_kineticEnergy)

add_executable(Test_leapfrog leapfrog.cpp)
target_link_libraries(Test_leapfrog ${Boost_LIBRARIES})
add_test(TestLeapfrogIntegration Test_leapfrog)

add_executable(Test_HMCProposal HMCProposal.cpp)
target_link_libraries(Test_HMCProposal ${Boost_LIBRARIES})
add_test(TestHamiltonianMonteCarloProposal Test_HMCProposal)
