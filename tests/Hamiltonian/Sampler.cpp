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

#define BOOST_TEST_MODULE Sampler
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>

BOOST_AUTO_TEST_CASE(sampler_finite_samples)
{
	typedef double RealType;
	typedef mcpack::utils::GaussPotentialEnergy<RealType> PotEngType;
	typedef PotEngType::RealVectorType RealVectorType;
	typedef PotEngType::RealMatrixType RealMatrixType;
	typedef RealMatrixType::Index IndexType;
	typedef mcpack::hamiltonian::GaussKineticEnergy<RealType> KinEngType;
	typedef mcpack::hamiltonian::LeapFrog<PotEngType,KinEngType> IntegratorType;
	typedef mcpack::hamiltonian::ClassicHMC<IntegratorType> HMCType;
	typedef mcpack::hamiltonian::IO_WriteAll<RealMatrixType> IOType;
	typedef mcpack::hamiltonian::RunCtrl_FiniteSamples<RealMatrixType> RCType;
	typedef mcpack::hamiltonian::Sampler<HMCType,IOType,RCType> SamplerType;

	
	//define the HMC
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

	//define the IO
	const std::string FileName("TestSampler.extract");
	IOType iowall(FileName);

	//define the Runtime Control
	const IndexType NumParas=10;
	const IndexType NSamples=1000;
	const IndexType NBurn=2000;
	const IndexType PacketSize=100;
	const std::string FileRoot("./TestSampler");
	const bool silent=false;

	RCType runctrl(NumParas,NSamples,PacketSize,NBurn,FileRoot,silent);

	SamplerType Smp(hmc,iowall,runctrl);

	Smp.Run();
}
