# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 2.8)
project(TestHamiltonianMonteCarloUtils CXX)

add_executable(Test_gaussLogPost gaussLogPost.cpp)
target_link_libraries(Test_gaussLogPost ${Boost_LIBRARIES})
add_test(TestGaussianLogPosterior Test_gaussLogPost)
