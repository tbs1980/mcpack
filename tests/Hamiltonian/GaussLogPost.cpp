/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

 
#define BOOST_TEST_MODULE GaussLogPost
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>



BOOST_AUTO_TEST_CASE(gaussian_log_posterior)
{
    typedef double RealType;
    typedef mcpack::utils::GaussPotentialEnergy<RealType> PotEngType;
    typedef PotEngType::RealVectorType RealVectorType;
    typedef PotEngType::RealMatrixType RealMatrixType;
    typedef RealMatrixType::Index IndexType;

    const IndexType N=100;

    RealVectorType mu=RealVectorType::Zero(N);
    RealMatrixType SigmaInv=RealMatrixType::Identity(N,N);
    RealVectorType q=RealVectorType::Random(N);

    PotEngType G(mu,SigmaInv);
    
    RealVectorType dq=SigmaInv*(mu-q);
    RealType val=-0.5*(mu-q).transpose()*SigmaInv*(mu-q);
    
    RealVectorType dqTest=RealVectorType::Zero(N);
    RealType valTest=0;
    G.Evaluate(q,valTest,dqTest);
    BOOST_CHECK_EQUAL(val,valTest);
    
    for(IndexType i=0;i<N;++i)
    {
        BOOST_CHECK_EQUAL(dq(i),dqTest(i));
    }
    
}

BOOST_AUTO_TEST_CASE(gaussian_log_posterior_diag_sigma_inv)
{
    typedef double RealType;
    typedef mcpack::utils::GaussPotentialEnergyDiag<RealType> PotEngType;
    typedef PotEngType::RealVectorType RealVectorType;
    typedef PotEngType::RealDiagMatrixType RealDiagMatrixType;
    typedef RealDiagMatrixType::Index IndexType;

    const IndexType N=100;

    RealVectorType mu=RealVectorType::Zero(N);
    RealVectorType q=RealVectorType::Random(N);
    RealDiagMatrixType SigmaInv(N);
    for(IndexType i=0;i<N;++i)
    {
        SigmaInv(i)=1;
    }

    PotEngType G(mu,SigmaInv);
    
    RealVectorType dq=SigmaInv.cwiseProduct(mu-q);
    RealType val=-0.5*(mu-q).transpose()*(SigmaInv.cwiseProduct(mu-q));

    RealVectorType dqTest=RealVectorType::Zero(N);
    RealType valTest=0;
    G.Evaluate(q,valTest,dqTest);

    
    BOOST_REQUIRE(fabs(val-valTest) < 1e-15);
    
    for(IndexType i=0;i<N;++i)
    {
        BOOST_REQUIRE( fabs(dq(i)-dqTest(i)) < 1e-15);
    }
}

