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

#define BOOST_TEST_MODULE ClassicHMC
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>


BOOST_AUTO_TEST_CASE(classic_hmc_10)
{
    typedef double RealType;
    typedef mcpack::utils::GaussPotentialEnergy<RealType> PotEngType;
    typedef PotEngType::RealVectorType RealVectorType;
    typedef PotEngType::RealMatrixType RealMatrixType;
    typedef RealMatrixType::Index IndexType;
    typedef mcpack::hamiltonian::GaussKineticEnergy<RealType> KinEngType;
    typedef mcpack::hamiltonian::LeapFrog<PotEngType,KinEngType> IntegratorType;
    typedef mcpack::hamiltonian::ClassicHMC<IntegratorType> HMCType;

    const IndexType N=10;

    RealVectorType mu=RealVectorType::Zero(N);
    RealMatrixType SigmaInv=RealMatrixType::Identity(N,N);
    RealMatrixType MInv=RealMatrixType::Identity(N,N);
    RealVectorType q0=RealVectorType::Random(N);

    const RealType eps=1;
    const IndexType Nsteps=10;

    PotEngType G(mu,SigmaInv);
    KinEngType K(MInv);

    IntegratorType Lp(G,K);

    HMCType hmc(Lp,eps,Nsteps,12346l,q0);

    const IndexType NSamples=1000;

    RealMatrixType Samples(NSamples,N);

    hmc.Generate(Samples);
    
    std::cout<<"Acceptace Rate= "<<hmc.GetAcceptanceRate()<<std::endl;

    BOOST_REQUIRE(0.91 < hmc.GetAcceptanceRate() and hmc.GetAcceptanceRate() < 0.95);

}


BOOST_AUTO_TEST_CASE(classic_hmc_100)
{
    typedef double RealType;
    typedef mcpack::utils::GaussPotentialEnergy<RealType> PotEngType;
    typedef PotEngType::RealVectorType RealVectorType;
    typedef PotEngType::RealMatrixType RealMatrixType;
    typedef RealMatrixType::Index IndexType;
    typedef mcpack::hamiltonian::GaussKineticEnergy<RealType> KinEngType;
    typedef mcpack::hamiltonian::LeapFrog<PotEngType,KinEngType> IntegratorType;
    typedef mcpack::hamiltonian::ClassicHMC<IntegratorType> HMCType;

    const IndexType N=100;

    RealVectorType mu=RealVectorType::Zero(N);
    RealMatrixType SigmaInv=RealMatrixType::Identity(N,N);
    RealMatrixType MInv=RealMatrixType::Identity(N,N);
    RealVectorType q0=RealVectorType::Random(N);

    const RealType eps=1;
    const IndexType Nsteps=10;

    PotEngType G(mu,SigmaInv);
    KinEngType K(MInv);

    IntegratorType Lp(G,K);

    HMCType hmc(Lp,eps,Nsteps,12346l,q0);

    const IndexType NSamples=1000;

    RealMatrixType Samples(NSamples,N);

    hmc.Generate(Samples);
    
    std::cout<<"Acceptace Rate= "<<hmc.GetAcceptanceRate()<<std::endl;

    BOOST_REQUIRE(0.79 < hmc.GetAcceptanceRate() and hmc.GetAcceptanceRate() < 0.83);

}

BOOST_AUTO_TEST_CASE(classic_hmc_diag_mass_mat)
{
    typedef double RealType;
    typedef mcpack::utils::GaussPotentialEnergy<RealType> PotEngType;
    typedef PotEngType::RealVectorType RealVectorType;
    typedef PotEngType::RealMatrixType RealMatrixType;
    typedef RealMatrixType::Index IndexType;
    typedef mcpack::hamiltonian::GaussKineticEnergyDiag<RealType> KinEngType;
    typedef KinEngType::RealDiagMatrixType RealDiagMatrixType;
    typedef mcpack::hamiltonian::LeapFrog<PotEngType,KinEngType> IntegratorType;
    typedef mcpack::hamiltonian::ClassicHMC<IntegratorType> HMCType;

    const IndexType N=100;

    RealVectorType mu=RealVectorType::Zero(N);
    RealMatrixType SigmaInv=RealMatrixType::Identity(N,N);
    RealVectorType q0=RealVectorType::Random(N);

    const RealType eps=1;
    const IndexType Nsteps=10;

    RealDiagMatrixType MInv(N);
    for(IndexType i=0;i<N;++i)
    {
        MInv(i)=1;
    }

    PotEngType G(mu,SigmaInv);
    KinEngType K(MInv);

    IntegratorType Lp(G,K);

    HMCType hmc(Lp,eps,Nsteps,12346l,q0);

    const IndexType NSamples=1000;

    RealMatrixType Samples(NSamples,N);

    hmc.Generate(Samples);
    
    std::cout<<"Acceptace Rate= "<<hmc.GetAcceptanceRate()<<std::endl;

    BOOST_REQUIRE(0.79 < hmc.GetAcceptanceRate() and hmc.GetAcceptanceRate() < 0.83);

}

BOOST_AUTO_TEST_CASE(classic_hmc_diag_mass_mat_sigma_inv)
{
    typedef double RealType;
    typedef mcpack::utils::GaussPotentialEnergyDiag<RealType> PotEngType;
    typedef PotEngType::RealVectorType RealVectorType;
    typedef PotEngType::RealDiagMatrixType RealDiagMatrixType;
    typedef RealDiagMatrixType::Index IndexType;
    typedef mcpack::hamiltonian::GaussKineticEnergyDiag<RealType> KinEngType;
    typedef KinEngType::RealDiagMatrixType RealDiagMatrixType;
    typedef mcpack::hamiltonian::LeapFrog<PotEngType,KinEngType> IntegratorType;
    typedef mcpack::hamiltonian::ClassicHMC<IntegratorType> HMCType;
    typedef Eigen::MatrixXd RealMatrixType;

    const IndexType N=1000;

    RealVectorType mu=RealVectorType::Zero(N);
    RealVectorType q0=RealVectorType::Random(N);

    const RealType eps=1;
    const IndexType Nsteps=10;

    RealDiagMatrixType MInv(N);
    RealDiagMatrixType SigmaInv(N);
   for(IndexType i=0;i<N;++i)
    {
        RealType s=2;
        SigmaInv(i)=s;
        MInv(i)=1./s;
    }

    PotEngType G(mu,SigmaInv);
    KinEngType K(MInv);

    IntegratorType Lp(G,K);

    HMCType hmc(Lp,eps,Nsteps,12346l,q0);

    const IndexType NSamples=1000;

    RealMatrixType Samples(NSamples,N);

    hmc.Generate(Samples);
    
    std::cout<<"Acceptace Rate= "<<hmc.GetAcceptanceRate()<<std::endl;

    BOOST_REQUIRE(0.56 < hmc.GetAcceptanceRate() and hmc.GetAcceptanceRate() < 0.58);

}

BOOST_AUTO_TEST_CASE(classic_hmc_2d_correlated)
{
    typedef double RealType;
    typedef mcpack::utils::GaussPotentialEnergy<RealType> PotEngType;
    typedef PotEngType::RealVectorType RealVectorType;
    typedef PotEngType::RealMatrixType RealMatrixType;
    typedef RealMatrixType::Index IndexType;
    typedef mcpack::hamiltonian::GaussKineticEnergy<RealType> KinEngType;
    typedef mcpack::hamiltonian::LeapFrog<PotEngType,KinEngType> IntegratorType;
    typedef mcpack::hamiltonian::ClassicHMC<IntegratorType> HMCType;

    const IndexType N=2;

    RealVectorType mu=RealVectorType::Zero(N);
    RealMatrixType SigmaInv=RealMatrixType::Identity(N,N);
    RealMatrixType MInv=RealMatrixType::Identity(N,N);
    RealVectorType q0=RealVectorType::Random(N);

    
    //assign the off diagonal terms
    SigmaInv(0,1) = 0.5;
    SigmaInv(1,0) = 0.5;
    //set MInv = (SigmaInv)^-1
    MInv = SigmaInv.inverse();
    

    std::cout<<"SigmaInv \n"<<SigmaInv<<"\n"<<std::endl;
    std::cout<<"MInv \n"<<MInv<<"\n"<<std::endl;

    
    const RealType eps=1;
    const IndexType Nsteps=10;

    PotEngType G(mu,SigmaInv);
    KinEngType K(MInv);

    IntegratorType Lp(G,K);

    HMCType hmc(Lp,eps,Nsteps,12346l,q0);

    const IndexType NSamples=1000;

    RealMatrixType Samples(NSamples,N);

    hmc.Generate(Samples);
    
    std::cout<<"Acceptace Rate= "<<hmc.GetAcceptanceRate()<<std::endl;

    BOOST_REQUIRE(0.98 < hmc.GetAcceptanceRate() and hmc.GetAcceptanceRate() < 0.99);
}

BOOST_AUTO_TEST_CASE(classic_hmc_10d_correlated)
{
    typedef double RealType;
    typedef mcpack::utils::GaussPotentialEnergy<RealType> PotEngType;
    typedef PotEngType::RealVectorType RealVectorType;
    typedef PotEngType::RealMatrixType RealMatrixType;
    typedef RealMatrixType::Index IndexType;
    typedef mcpack::hamiltonian::GaussKineticEnergy<RealType> KinEngType;
    typedef mcpack::hamiltonian::LeapFrog<PotEngType,KinEngType> IntegratorType;
    typedef mcpack::hamiltonian::ClassicHMC<IntegratorType> HMCType;

    const IndexType N=10;

    RealVectorType mu=RealVectorType::Zero(N);
    RealVectorType q0=RealVectorType::Random(N);

    RealMatrixType Mat=RealMatrixType::Random(N,N);
    RealMatrixType Sigma=Mat*Mat.transpose();;
    RealMatrixType SigmaInv=Sigma.inverse();
    RealMatrixType MInv=Sigma;

    std::cout<<"SigmaInv \n"<<SigmaInv<<"\n"<<std::endl;
    std::cout<<"MInv \n"<<MInv<<"\n"<<std::endl;

    
    const RealType eps=1;
    const IndexType Nsteps=10;

    PotEngType G(mu,SigmaInv);
    KinEngType K(MInv);

    IntegratorType Lp(G,K);

    HMCType hmc(Lp,eps,Nsteps,12346l,q0);

    const IndexType NSamples=1000;

    RealMatrixType Samples(NSamples,N);

    hmc.Generate(Samples);
    
    std::cout<<"Acceptace Rate= "<<hmc.GetAcceptanceRate()<<std::endl;

    BOOST_REQUIRE(0.93 < hmc.GetAcceptanceRate() and hmc.GetAcceptanceRate() < 0.94);
}
