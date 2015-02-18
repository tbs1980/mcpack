/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#define BOOST_TEST_MODULE RuntimeControl
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>

template<typename realScalarType>
void test_finite_samples(void)
{
    typedef Eigen::Matrix<realScalarType, Eigen::Dynamic, Eigen::Dynamic> realMatrixType;
    typedef typename realMatrixType::Index IndexType;
    typedef mcpack::hamiltonian::runCtrlFiniteSamples<realMatrixType> rCType;
    typedef mcpack::utils::randomSTD<realScalarType> randVarGenType;

    randVarGenType rngVar(0);
    std::stringstream rngState;
    rngVar.getState(rngState);
    realScalarType accRate=0.9;

    const IndexType numParas=10;
    const IndexType numsamples=1000;
    const IndexType numBurn=200;
    const IndexType packetSize=100;
    const std::string fileRoot("./testRunCtrl");
    const bool silent=false;
    rCType runctrl(numParas,numsamples,packetSize,numBurn,fileRoot,silent);

    BOOST_REQUIRE(runctrl.continueSampling());

    for(IndexType i=0;i<(numsamples+numBurn)/packetSize-1;++i)
    {
        realMatrixType samples=realMatrixType::Random(packetSize,numParas);
        runctrl.save(samples,rngState,accRate);
        BOOST_REQUIRE(runctrl.continueSampling());
    }

    realMatrixType samples=realMatrixType::Random(packetSize,numParas);
    runctrl.save(samples,rngState,accRate);
    BOOST_REQUIRE( !runctrl.continueSampling() );

}

BOOST_AUTO_TEST_CASE(finite_samples)
{
    test_finite_samples<float>();
    test_finite_samples<double>();
    test_finite_samples<long double>();
}
