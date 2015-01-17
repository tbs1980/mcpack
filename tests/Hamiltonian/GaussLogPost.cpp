/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#define BOOST_TEST_MODULE GaussLogPost
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>

template<typename realType>
void test_gaussian_log_posterior(void)
{
    typedef typename mcpack::utils::gaussPotentialEnergy<realType> potEngType;
    typedef typename potEngType::realVectorType realVectorType;
    typedef typename potEngType::realMatrixType realMatrixType;
    typedef typename realMatrixType::Index indexType;

    const indexType N=100;

    realVectorType mu=realVectorType::Zero(N);
    realMatrixType sigmaInv=realMatrixType::Identity(N,N);
    realVectorType q=realVectorType::Random(N);

    potEngType G(mu,sigmaInv);

    realVectorType dq=sigmaInv*(mu-q);
    realType val=-0.5*(mu-q).transpose()*sigmaInv*(mu-q);

    realVectorType dqTest=realVectorType::Zero(N);
    realType valTest=0;
    G.evaluate(q,valTest,dqTest);
    BOOST_CHECK_EQUAL(val,valTest);

    for(indexType i=0;i<N;++i)
    {
        BOOST_CHECK_EQUAL(dq(i),dqTest(i));
    }
}

template<typename realType>
void test_gaussian_log_posterior_diag_sigma_inv(void)
{
    typedef typename mcpack::utils::gaussPotentialEnergyDiag<realType> potEngType;
    typedef typename potEngType::realVectorType realVectorType;
    typedef typename potEngType::realDiagMatrixType realDiagMatrixType;
    typedef typename realDiagMatrixType::Index indexType;

    const indexType N=100;

    realVectorType mu=realVectorType::Zero(N);
    realVectorType q=realVectorType::Random(N);
    realDiagMatrixType sigmaInv(N);
    for(indexType i=0;i<N;++i)
    {
        sigmaInv(i)=1;
    }

    potEngType G(mu,sigmaInv);

    realVectorType dq=sigmaInv.cwiseProduct(mu-q);
    realType val=-0.5*(mu-q).transpose()*(sigmaInv.cwiseProduct(mu-q));

    realVectorType dqTest=realVectorType::Zero(N);
    realType valTest=0;
    G.evaluate(q,valTest,dqTest);


    BOOST_REQUIRE(fabs(val-valTest) < 1e-15);

    for(indexType i=0;i<N;++i)
    {
        BOOST_REQUIRE( fabs(dq(i)-dqTest(i)) < 1e-15);
    }
}

BOOST_AUTO_TEST_CASE(gaussian_log_posterior)
{
    test_gaussian_log_posterior<float>();
    test_gaussian_log_posterior<double>();
    test_gaussian_log_posterior<long double>();
}

BOOST_AUTO_TEST_CASE(gaussian_log_posterior_diag_sigma_inv)
{
    test_gaussian_log_posterior_diag_sigma_inv<float>();
    test_gaussian_log_posterior_diag_sigma_inv<double>();
    test_gaussian_log_posterior_diag_sigma_inv<long double>();
}
