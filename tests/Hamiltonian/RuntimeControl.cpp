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

#define BOOST_TEST_MODULE RuntimeControl
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>

BOOST_AUTO_TEST_CASE(finite_samples)
{
	typedef double RealType;
	typedef Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic> RealMatrixType;
	typedef RealMatrixType::Index IndexType;
	typedef mcpack::hamiltonian::RunCtrl_FiniteSamples<RealMatrixType> RCType;

	const IndexType NParas=10;
	const IndexType NSamples=1000;
	const IndexType NBurn=200;
	const IndexType PacketSize=100;
	RCType runctrl(NParas,NSamples,PacketSize,NBurn);

	BOOST_REQUIRE(runctrl.Continue());

	for(IndexType i=0;i<(NSamples+NBurn)/PacketSize-1;++i)
	{
		RealMatrixType Samples=RealMatrixType::Random(PacketSize,NParas);
		runctrl.Add(Samples);
		BOOST_REQUIRE(runctrl.Continue());
	}

	RealMatrixType Samples=RealMatrixType::Random(PacketSize,NParas);
	runctrl.Add(Samples);
	BOOST_REQUIRE( !runctrl.Continue() );

}