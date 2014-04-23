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

#define BOOST_TEST_MODULE ClassicHMC
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>

BOOST_AUTO_TEST_CASE(classic_hmc_10)
{
	typedef double RealType;
	typedef mcpack::utils::GaussPotentialEnergy<RealType> PotEngType;
	typedef PotEngType::RealVectorType RealVectorType;
	typedef PotEngType::RealMatrixType RealMatrixType;
	typedef RealMatrixType::Index IndexType;
	typedef mcpack::hamiltonian::GaussKineticEnergy<RealType> KinEngType;
	typedef mcpack::hamiltonian::LeapFrog<PotEngType,KinEngType> IntegratorType;
	typedef mcpack::hamiltonian::ClassicHMC<IntegratorType> HMCType;

	const IndexType N=10;

	RealVectorType mu=RealVectorType::Zero(N);
	RealMatrixType SigmaInv=RealMatrixType::Identity(N,N);
	RealMatrixType MInv=RealMatrixType::Identity(N,N);
	RealVectorType q0=RealVectorType::Random(N);

	const RealType eps=1;
	const IndexType Nsteps=10;

	PotEngType G(mu,SigmaInv);
	KinEngType K(MInv);

	IntegratorType Lp(G,K);

	HMCType hmc(Lp,eps,Nsteps,12346l,q0);

	const IndexType NSamples=1000;

	RealMatrixType Samples(NSamples,N);

	hmc.Generate(Samples);
	
	std::cout<<"Acceptace Rate= "<<hmc.GetAcceptanceRate()<<std::endl;

	BOOST_REQUIRE(0.91 < hmc.GetAcceptanceRate() and hmc.GetAcceptanceRate() < 0.95);

}



BOOST_AUTO_TEST_CASE(classic_hmc_100)
{
	typedef double RealType;
	typedef mcpack::utils::GaussPotentialEnergy<RealType> PotEngType;
	typedef PotEngType::RealVectorType RealVectorType;
	typedef PotEngType::RealMatrixType RealMatrixType;
	typedef RealMatrixType::Index IndexType;
	typedef mcpack::hamiltonian::GaussKineticEnergy<RealType> KinEngType;
	typedef mcpack::hamiltonian::LeapFrog<PotEngType,KinEngType> IntegratorType;
	typedef mcpack::hamiltonian::ClassicHMC<IntegratorType> HMCType;

	const IndexType N=100;

	RealVectorType mu=RealVectorType::Zero(N);
	RealMatrixType SigmaInv=RealMatrixType::Identity(N,N);
	RealMatrixType MInv=RealMatrixType::Identity(N,N);
	RealVectorType q0=RealVectorType::Random(N);

	const RealType eps=1;
	const IndexType Nsteps=10;

	PotEngType G(mu,SigmaInv);
	KinEngType K(MInv);

	IntegratorType Lp(G,K);

	HMCType hmc(Lp,eps,Nsteps,12346l,q0);

	const IndexType NSamples=1000;

	RealMatrixType Samples(NSamples,N);

	hmc.Generate(Samples);
	
	std::cout<<"Acceptace Rate= "<<hmc.GetAcceptanceRate()<<std::endl;

	BOOST_REQUIRE(0.79 < hmc.GetAcceptanceRate() and hmc.GetAcceptanceRate() < 0.83);

}