/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
 
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
    typedef typename mcpack::utils::RandomVariateGenerator<RealType> RandVarGenType;
    typedef mcpack::hamiltonian::HMCProposal<IntegratorType,RandVarGenType> HMCProposalType;

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

    RealType dH=0;
    RandVarGenType RVGen(12346l);

    prop.Propose(q0,dH,RVGen);
}
