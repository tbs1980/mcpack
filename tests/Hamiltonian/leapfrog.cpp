/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#define BOOST_TEST_MODULE leapfrog
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>

template<typename realScalarType>
void test_Identiy_Matrix_Gauss(void)
{
    typedef mcpack::utils::gaussPotentialEnergy<realScalarType> potEngType;
    typedef typename potEngType::realVectorType realVectorType;
    typedef typename potEngType::realMatrixType realMatrixType;
    typedef typename realMatrixType::Index indexType;
    typedef mcpack::hamiltonian::gaussKineticEnergy<realScalarType> kinEngType;
    typedef mcpack::hamiltonian::leapfrog<potEngType,kinEngType> leapfrogIntegratorType;

    const indexType N=100;

    realVectorType q=realVectorType::Random(N);
    realVectorType dq=realVectorType::Zero(N);
    realVectorType p=realVectorType::Random(N);

    realVectorType mu=realVectorType::Zero(N);
    realMatrixType sigmaInv=realMatrixType::Identity(N,N);
    realMatrixType mInv=realMatrixType::Identity(N,N);

    const realScalarType eps=1;
    const indexType nSteps=10;

    // after one iteration we should get
    // q(t+e) = q(t) + e*p(t) - 0.5*e^2*q(t)
    // p(t+e) = (1-0.5*e^2)*p(t)+ (0.25*e^3-e)*q(t)

    realVectorType qtest(q);
    realVectorType ptest(p);

    for(indexType i=0;i<nSteps;++i)
    {
        realVectorType qtemp = qtest + eps*ptest - 0.5*eps*eps*qtest;
        realVectorType ptemp = (1-0.5*eps*eps)*ptest + (0.25*eps*eps*eps-eps)*qtest;
        qtest=qtemp;
        ptest=ptemp;
    }


    potEngType G(mu,sigmaInv);
    kinEngType K(mInv);

    leapfrogIntegratorType lP(G,K);

    realScalarType dH=0;

    lP.integrate(q,p,eps,nSteps,dH);

    realScalarType meps=std::numeric_limits<realScalarType>::epsilon();

    for(indexType i=0;i<N;++i)
    {
        BOOST_REQUIRE( std::abs(q(i)-qtest(i)) < 1e2*meps );
        BOOST_REQUIRE( std::abs(p(i)-ptest(i)) < 1e2*meps );
    }

}

BOOST_AUTO_TEST_CASE(Identiy_Matrix_Gauss)
{
    test_Identiy_Matrix_Gauss<float>();
    test_Identiy_Matrix_Gauss<double>();
    test_Identiy_Matrix_Gauss<long double>();
}
