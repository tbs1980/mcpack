/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#define BOOST_TEST_MODULE KineticEnergy
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>

BOOST_AUTO_TEST_CASE(Nd_General_Gauss_Test)
{
    typedef double RealType;
    typedef mcpack::hamiltonian::GaussKineticEnergy<RealType> KinEngType;
    typedef KinEngType::RealVectorType RealVectorType;
    typedef KinEngType::RealMatrixType RealMatrixType;
    typedef RealMatrixType::Index IndexType;

    const IndexType N=100;

    RealVectorType x=RealVectorType::Random(N);
    RealVectorType g=RealVectorType::Zero(N);
    RealMatrixType MInv=RealMatrixType::Identity(N,N);
    RealType val=0;
    KinEngType phi(MInv);

    phi.Evaluate(x,val,g);

    RealType valTest=-0.5*x.transpose()*MInv*x;

    BOOST_CHECK_EQUAL(val,valTest);

    for(IndexType i=0;i<N;++i)
    {
        BOOST_CHECK_EQUAL(g(i),-x(i));
    }

    BOOST_CHECK_EQUAL(N,phi.NDim());
}

BOOST_AUTO_TEST_CASE(Nd_Diag_Gauss_Test)
{
    typedef double RealType;
    typedef mcpack::hamiltonian::GaussKineticEnergyDiag<RealType> KinEngType;
    typedef KinEngType::RealVectorType RealVectorType;
    typedef KinEngType::RealDiagMatrixType RealDiagMatrixType;
    typedef RealDiagMatrixType::Index IndexType;

    const IndexType N=1000000;

    RealVectorType x=RealVectorType::Random(N);
    RealVectorType g=RealVectorType::Zero(N);
    RealDiagMatrixType MInv(N);
    for(IndexType i=0;i<N;++i)
    {
        MInv(i)=1;
    }
    RealType val=0;
    KinEngType phi(MInv);

    phi.Evaluate(x,val,g);

    RealType valTest=-0.5*x.transpose()*(MInv.cwiseProduct(x));

    BOOST_CHECK_EQUAL(val,valTest);

    for(IndexType i=0;i<N;++i)
    {
        BOOST_CHECK_EQUAL(g(i),-x(i));
    }

    BOOST_CHECK_EQUAL(N,phi.NDim());
}