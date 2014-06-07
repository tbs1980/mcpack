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

