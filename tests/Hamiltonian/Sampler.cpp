/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#define BOOST_TEST_MODULE Sampler
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>

BOOST_AUTO_TEST_CASE(sampler_finite_samples)
{
    typedef double RealType;
    typedef mcpack::utils::GaussPotentialEnergy<RealType> PotEngType;
    typedef PotEngType::RealVectorType RealVectorType;
    typedef PotEngType::RealMatrixType RealMatrixType;
    typedef RealMatrixType::Index IndexType;
    typedef mcpack::hamiltonian::GaussKineticEnergy<RealType> KinEngType;
    typedef mcpack::hamiltonian::LeapFrog<PotEngType,KinEngType> IntegratorType;
    typedef mcpack::hamiltonian::ClassicHMC<IntegratorType> HMCType;
    typedef mcpack::hamiltonian::IO_WriteAll<RealMatrixType> IOType;
    typedef mcpack::hamiltonian::RunCtrl_FiniteSamples<RealMatrixType> RCType;
    typedef mcpack::hamiltonian::Sampler<HMCType,IOType,RCType> SamplerType;

    
    //define the HMC
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

    //define the IO
    const std::string FileName("TestSampler.extract");
    IOType iowall(FileName);

    //define the Runtime Control
    const IndexType NumParas=10;
    const IndexType NSamples=1000;
    const IndexType NBurn=2000;
    const IndexType PacketSize=100;
    const std::string FileRoot("./TestSampler");
    const bool silent=false;

    RCType runctrl(NumParas,NSamples,PacketSize,NBurn,FileRoot,silent);

    SamplerType Smp(hmc,iowall,runctrl);

    Smp.Run();
}

BOOST_AUTO_TEST_CASE(sampler_finite_samples_diag_KE_and_PE)
{
    typedef double RealType;
    typedef mcpack::utils::GaussPotentialEnergyDiag<RealType> PotEngType;
    typedef PotEngType::RealVectorType RealVectorType;
    typedef PotEngType::RealDiagMatrixType RealDiagMatrixType;
    typedef RealDiagMatrixType::Index IndexType;
    typedef Eigen::MatrixXd RealMatrixType;
    typedef mcpack::hamiltonian::GaussKineticEnergyDiag<RealType> KinEngType;
    typedef mcpack::hamiltonian::LeapFrog<PotEngType,KinEngType> IntegratorType;
    typedef mcpack::hamiltonian::ClassicHMC<IntegratorType> HMCType;
    typedef mcpack::hamiltonian::IO_WriteAll<RealMatrixType> IOType;
    typedef mcpack::hamiltonian::RunCtrl_FiniteSamples<RealMatrixType> RCType;
    typedef mcpack::hamiltonian::Sampler<HMCType,IOType,RCType> SamplerType;

    
    //define the HMC
    const IndexType N=10;
    RealVectorType mu=RealVectorType::Zero(N);
    RealDiagMatrixType SigmaInv(N);
    RealDiagMatrixType MInv(N);
    RealVectorType q0=RealVectorType::Random(N);
    for(IndexType i=0;i<N;++i)
    {
        SigmaInv(i)=1;
        MInv(i)=1;
    }

    const RealType eps=1;
    const IndexType Nsteps=10;

    PotEngType G(mu,SigmaInv);
    KinEngType K(MInv);

    IntegratorType Lp(G,K);

    HMCType hmc(Lp,eps,Nsteps,12346l,q0);

    //define the IO
    const std::string FileName("TestSamplerDiag.extract");
    IOType iowall(FileName);

    //define the Runtime Control
    const IndexType NumParas=10;
    const IndexType NSamples=1000;
    const IndexType NBurn=2000;
    const IndexType PacketSize=100;
    const std::string FileRoot("./TestSamplerDiag");
    const bool silent=false;

    RCType runctrl(NumParas,NSamples,PacketSize,NBurn,FileRoot,silent);

    SamplerType Smp(hmc,iowall,runctrl);

    Smp.Run();
}



