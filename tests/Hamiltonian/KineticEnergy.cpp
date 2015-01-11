/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#define BOOST_TEST_MODULE KineticEnergy
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <mcpack/CoreHeaders.hpp>

BOOST_AUTO_TEST_CASE(Nd_General_Gauss_Test)
{
    typedef double realType;
    typedef mcpack::hamiltonian::gaussKineticEnergy<realType> kinEngType;
    typedef kinEngType::realVectorType realVectorType;
    typedef kinEngType::realMatrixType realMatrixType;
    typedef realMatrixType::Index indexType;

    const indexType N=100;

    realVectorType x=realVectorType::Random(N);
    realVectorType g=realVectorType::Zero(N);
    realMatrixType mInv=realMatrixType::Identity(N,N);
    realType val=0;
    kinEngType phi(mInv);

    phi.evaluate(x,val,g);

    realType valTest=-0.5*x.transpose()*mInv*x;

    BOOST_CHECK_EQUAL(val,valTest);

    for(indexType i=0;i<N;++i)
    {
        BOOST_CHECK_EQUAL(g(i),-x(i));
    }

    BOOST_CHECK_EQUAL(N,phi.numDims());
}


BOOST_AUTO_TEST_CASE(Nd_Diag_Gauss_Test)
{
    typedef double realType;
    typedef mcpack::hamiltonian::gaussKineticEnergyDiag<realType> kinEngType;
    typedef kinEngType::realVectorType realVectorType;
    typedef kinEngType::realDiagMatrixType realDiagMatrixType;
    typedef realDiagMatrixType::Index indexType;

    const indexType N=1000000;

    realVectorType x=realVectorType::Random(N);
    realVectorType g=realVectorType::Zero(N);
    realDiagMatrixType mInv(N);
    for(indexType i=0;i<N;++i)
    {
        mInv(i)=1;
    }
    realType val=0;
    kinEngType phi(mInv);

    phi.evaluate(x,val,g);

    realType valTest=-0.5*x.transpose()*(mInv.cwiseProduct(x));

    BOOST_CHECK_EQUAL(val,valTest);

    for(indexType i=0;i<N;++i)
    {
        BOOST_CHECK_EQUAL(g(i),-x(i));
    }

    BOOST_CHECK_EQUAL(N,phi.numDims());
}
