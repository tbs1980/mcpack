/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#define BOOST_TEST_MODULE IO
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>

template<typename realScalarType>
void test_finite_samples(void)
{
    typedef Eigen::Matrix<realScalarType, Eigen::Dynamic, Eigen::Dynamic> realMatrixType;
    typedef typename realMatrixType::Index indexType;
    typedef mcpack::hamiltonian::IOWriteAllParams<realMatrixType> IOType;

    const indexType numParams=10;
    const indexType packetSize=100;
    const std::string fileName("testIO.dat");
    //clear the contents of TestIO.dat if previous versions exist
    std::ofstream of;
    of.open(fileName.c_str(),std::ios::trunc);
    of.close();

    IOType iowall(fileName);


    realMatrixType Samples=realMatrixType::Random(packetSize,numParams);

    iowall.write(Samples);

/*
    realMatrixType SamplesTest=realMatrixType::Zero(packetSize,numParams);

    const std::string separation(",");
    mcpack::utils::ReadMatrixFromTextFile(SamplesTest,fileName,separation);

    for(indexType i=0;i<packetSize;++i)
    {
        for(indexType j=0;j<numParams;++j)
        {
            realScalarType rel_diff=( Samples(i,j)-SamplesTest(i,j) )/Samples(i,j);
            BOOST_REQUIRE(fabs( rel_diff )  < 1e-5);
        }
    }
*/
}

BOOST_AUTO_TEST_CASE(finite_samples)
{
    test_finite_samples<float>();
    test_finite_samples<double>();
    test_finite_samples<long double>();
}
