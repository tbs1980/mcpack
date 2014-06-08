/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
    //clear the contents of TestIO.dat if previous versions exist
    std::ofstream of;
    of.open(FileName.c_str(),std::ios::trunc);
    of.close();

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