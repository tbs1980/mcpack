/* 
* 
* Copyright (C) 2014 Sreekumar Thaithara Balan
* 
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or (at
* your option) any later version.
* 
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

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

#include "Source/Utils/Macros.hpp"
#include "Source/Utils/Random.hpp"
#include "Source/Utils/GaussLogPost.hpp"
#include "Source/Utils/TextDataIO.hpp"

#include "Source/Hamiltonian/KineticEnergy.hpp"
#include "Source/Hamiltonian/Integrator.hpp"
#include "Source/Hamiltonian/ClassicHMC.hpp"
#include "Source/Hamiltonian/IO.hpp"
#include "Source/Hamiltonian/RuntimeControl.hpp"
#include "Source/Hamiltonian/Sampler.hpp"

#endif //MCPACK_COREHEADERS_HPP
