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

#define BOOST_TEST_MODULE Integration
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>

BOOST_AUTO_TEST_CASE(Identiy_Matrix_Gauss)
{
	typedef double RealType;
	typedef mcpack::utils::GaussPotentialEnergy<RealType> PotEngType;
	typedef PotEngType::RealVectorType RealVectorType;
	typedef PotEngType::RealMatrixType RealMatrixType;
	typedef RealMatrixType::Index IndexType;
	typedef mcpack::hamiltonian::GaussKineticEnergy<RealType> KinEngType;
	typedef mcpack::hamiltonian::LeapFrog<PotEngType,KinEngType> IntegratorType;

	const IndexType N=100;

	RealVectorType q=RealVectorType::Random(N);
	RealVectorType dq=RealVectorType::Zero(N);
	RealVectorType p=RealVectorType::Random(N);

	RealVectorType mu=RealVectorType::Zero(N);
	RealMatrixType SigmaInv=RealMatrixType::Identity(N,N);
	RealMatrixType MInv=RealMatrixType::Identity(N,N);

	const RealType eps=1;
	const IndexType Nsteps=10;

	// after one iteration we should get 
	// q(t+e) = q(t) + e*p(t) - 0.5*e^2*q(t)
	// p(t+e) = (1-0.5*e^2)*p(t)+ (0.25*e^3-e)*q(t)

	RealVectorType qtest(q);
	RealVectorType ptest(p);

	for(IndexType i=0;i<Nsteps;++i)
	{
		RealVectorType qtemp = qtest + eps*ptest - 0.5*eps*eps*qtest;
		RealVectorType ptemp = (1-0.5*eps*eps)*ptest + (0.25*eps*eps*eps-eps)*qtest;
		qtest=qtemp;
		ptest=ptemp;
	}


	PotEngType G(mu,SigmaInv);
	KinEngType K(MInv);

	IntegratorType Lp(G,K);

	RealType dH=0;

	Lp.Integrate(q,p,eps,Nsteps,dH);

	RealType meps=std::numeric_limits<RealType>::epsilon();

	for(IndexType i=0;i<N;++i)
	{
		BOOST_REQUIRE(fabs(q(i)-qtest(i)) < 1e2*meps );
		BOOST_REQUIRE(fabs(p(i)-ptest(i)) < 1e2*meps );
	}

}
