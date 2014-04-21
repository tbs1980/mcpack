cmake_minimum_required(VERSION 2.8)
project(MCPACKUnitTests CXX)
if(UNIT_TESTS)
	enable_testing()
endif()

add_subdirectory(Hamiltonian)
