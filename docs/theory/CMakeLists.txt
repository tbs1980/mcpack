# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

cmake_minimum_required(VERSION 2.8)
project(MCPACK_Theory NONE)
include(${CMAKE_SOURCE_DIR}/cmake/UseLATEX.cmake)
add_latex_document(mcpack.tex
	INPUTS include/hamiltonian.tex
	BIBFILES mcpack.bib 
#	IMAGE_DIRS images 
	DEFAULT_PDF)
