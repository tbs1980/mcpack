/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
    typedef mcpack::utils::RandomVariateGenerator<RealType> RandVarGenType;

    RandVarGenType RngVar(0);
    std::stringstream RngState;
    RngVar.GetState(RngState);
    RealType AccRate=0.9;

    const IndexType NParas=10;
    const IndexType NSamples=1000;
    const IndexType NBurn=200;
    const IndexType PacketSize=100;
    const std::string FileRoot("./TestRunCtrl");
    const bool silent=false;
    RCType runctrl(NParas,NSamples,PacketSize,NBurn,FileRoot,silent);

    BOOST_REQUIRE(runctrl.Continue());

    for(IndexType i=0;i<(NSamples+NBurn)/PacketSize-1;++i)
    {
        RealMatrixType Samples=RealMatrixType::Random(PacketSize,NParas);
        runctrl.Save(Samples,RngState,AccRate);
        BOOST_REQUIRE(runctrl.Continue());
    }

    RealMatrixType Samples=RealMatrixType::Random(PacketSize,NParas);
    runctrl.Save(Samples,RngState,AccRate);
    BOOST_REQUIRE( !runctrl.Continue() );

}