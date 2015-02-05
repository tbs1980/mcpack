/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#define BOOST_TEST_MODULE KineticEnergy
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <mcpack/CoreHeaders.hpp>

template<typename realScalarType>
void test_Nd_General_Gauss(void)
{
    typedef typename mcpack::hamiltonian::gaussKineticEnergy<realScalarType> kinEngType;
    typedef typename kinEngType::realVectorType realVectorType;
    typedef typename kinEngType::realMatrixType realMatrixType;
    typedef typename realMatrixType::Index indexType;

    const indexType N=100;

    realVectorType x=realVectorType::Random(N);
    realVectorType g=realVectorType::Zero(N);
    realMatrixType mInv=realMatrixType::Identity(N,N);
    realScalarType val=0;
    kinEngType phi(mInv);

    phi.evaluate(x,val,g);

    realScalarType valTest=-0.5*x.transpose()*mInv*x;

    BOOST_CHECK_EQUAL(val,valTest);

    for(indexType i=0;i<N;++i)
    {
        BOOST_CHECK_EQUAL(g(i),-x(i));
    }

    BOOST_CHECK_EQUAL(N,phi.numDims());
}

template<typename realScalarType>
void test_Nd_Diag_Gauss(void)
{
    typedef typename mcpack::hamiltonian::gaussKineticEnergyDiag<realScalarType> kinEngType;
    typedef typename kinEngType::realVectorType realVectorType;
    typedef typename kinEngType::realDiagMatrixType realDiagMatrixType;
    typedef typename realDiagMatrixType::Index indexType;

    const indexType N=1000000;

    realVectorType x=realVectorType::Random(N);
    realVectorType g=realVectorType::Zero(N);
    realDiagMatrixType mInv(N);
    for(indexType i=0;i<N;++i)
    {
        mInv(i)=1;
    }
    realScalarType val=0;
    kinEngType phi(mInv);

    phi.evaluate(x,val,g);

    realScalarType valTest=-0.5*x.transpose()*(mInv.cwiseProduct(x));

    BOOST_CHECK_EQUAL(val,valTest);

    for(indexType i=0;i<N;++i)
    {
        BOOST_CHECK_EQUAL(g(i),-x(i));
    }

    BOOST_CHECK_EQUAL(N,phi.numDims());
}

BOOST_AUTO_TEST_CASE(Nd_General_Gauss_Test)
{
    test_Nd_General_Gauss<float>();
    test_Nd_General_Gauss<double>();
    test_Nd_General_Gauss<long double>();
}


BOOST_AUTO_TEST_CASE(Nd_Diag_Gauss_Test)
{
    test_Nd_Diag_Gauss<float>();
    test_Nd_Diag_Gauss<double>();
    test_Nd_Diag_Gauss<long double>();
}
