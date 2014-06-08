/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <mcpack/MPICoreHeaders.hpp>

int main()
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
    typedef mcpack::hamiltonian::Mpi_Sampler<HMCType,IOType,RCType> MPISamplerType;

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
    const std::string FileName("TestMPISampler.extract");
    IOType iowall(FileName);

    //define the Runtime Control
    const IndexType NumParas=10;
    const IndexType NSamples=1000;
    const IndexType NBurn=200;
    const IndexType PacketSize=100;
    const std::string FileRoot("./TestMPISampler");
    const bool silent=false;
    RCType runctrl(NumParas,NSamples,PacketSize,NBurn,FileRoot,silent);

    MPISamplerType MPISmp(hmc,iowall,runctrl);

    MPISmp.Run();

    return 0;
}
