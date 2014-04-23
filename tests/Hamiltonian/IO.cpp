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

#define BOOST_TEST_MODULE IO
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>

BOOST_AUTO_TEST_CASE(finite_samples)
{
	typedef double RealType;
	typedef Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic> RealMatrixType;
	typedef RealMatrixType::Index IndexType;
	typedef mcpack::hamiltonian::IO_WriteAll<RealMatrixType> IOType;

	const IndexType NParas=10;
	const IndexType PacketSize=100;
	const std::string FileName("TestIO.dat");
	IOType iowall(FileName);

	RealMatrixType Samples=RealMatrixType::Random(PacketSize,NParas);

	iowall.Write(Samples);

	RealMatrixType SamplesTest=RealMatrixType::Zero(PacketSize,NParas);

	const std::string separation(",");
	mcpack::utils::ReadMatrixFromTextFile(SamplesTest,FileName,separation);

	for(IndexType i=0;i<PacketSize;++i)
	{
		for(IndexType j=0;j<NParas;++j)
		{
			RealType rel_diff=( Samples(i,j)-SamplesTest(i,j) )/Samples(i,j);
			BOOST_REQUIRE(fabs( rel_diff )  < 1e-5);
		}
	}

}