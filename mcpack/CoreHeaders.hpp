/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_COREHEADERS_HPP
#define MCPACK_COREHEADERS_HPP

#include <cstdlib>
#include <cmath>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <type_traits>
#include <type_traits>
#include <typeinfo>
#include <random>
#include <algorithm>
#include <exception>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/assert.hpp>

#include "Source/utils/randomSTD.hpp"
#include "Source/utils/gaussLogPost.hpp"

/*
#include "Source/Utils/Macros.hpp"
#include "Source/Utils/TextDataIO.hpp"
*/

#include "Source/Hamiltonian/kineticEnergy.hpp"
#include "Source/Hamiltonian/leapfrog.hpp"
#include "Source/Hamiltonian/HMCProposal.hpp"
#include "Source/Hamiltonian/classicHMC.hpp"
#include "Source/Hamiltonian/runtimeControl.hpp"

/*

#include "Source/Hamiltonian/IO.hpp"
#include "Source/Hamiltonian/Sampler.hpp"
*/

#endif //MCPACK_COREHEADERS_HPP
