/* 
* 
* Copyright (C) 2014 Sreekumar Thaithara Balan
* 
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or (at
* your option) any later version.
* 
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

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