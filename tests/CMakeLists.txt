# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 2.8)
project(MCPACKUnitTests CXX)
if(UNIT_TESTS)
    enable_testing()
endif()

add_subdirectory(Hamiltonian)
