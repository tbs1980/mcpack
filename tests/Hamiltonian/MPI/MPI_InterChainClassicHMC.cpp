/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <mcpack/MPICoreHeaders.hpp>

int main(void)
{
    typedef double RealType;
    typedef mcpack::utils::GaussPotentialEnergy<RealType> PotEngType;
    typedef PotEngType::RealVectorType RealVectorType;
    typedef PotEngType::RealMatrixType RealMatrixType;
    typedef RealMatrixType::Index IndexType;
    typedef mcpack::hamiltonian::GaussKineticEnergy<RealType> KinEngType;
    typedef mcpack::hamiltonian::LeapFrog<PotEngType,KinEngType> IntegratorType;
    typedef mcpack::utils::RandomVariateGenerator<RealType> RandVarGenType;
    typedef mcpack::hamiltonian::HMCProposal<IntegratorType,RandVarGenType> HMCProposalType;
    //typedef mcpack::hamiltonian::ClassicHMC<HMCProposalType> HMCType;
    typedef mcpack::hamiltonian::MPIInterChainClassicHMC<HMCProposalType> HMCType;

    
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

    HMCProposalType prop(Lp,eps,Nsteps);

    HMCType hmc(prop,q0,12346l);

    
    const IndexType NSamples=1000;

    RealMatrixType Samples(NSamples,N);

    hmc.Generate(Samples);
    
    std::cout<<"Acceptace Rate= "<<hmc.GetAcceptanceRate()<<std::endl;

	return 0;
}