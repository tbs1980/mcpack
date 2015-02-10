/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#define BOOST_TEST_MODULE ClassicHMC
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <mcpack/CoreHeaders.hpp>

BOOST_AUTO_TEST_CASE(classic_hmc_10)
{
    typedef double realScalarType;
    typedef mcpack::utils::gaussPotentialEnergy<realScalarType> potEngType;
    typedef potEngType::realVectorType realVectorType;
    typedef potEngType::realMatrixType realMatrixType;
    typedef realMatrixType::Index indexType;
    typedef mcpack::hamiltonian::gaussKineticEnergy<realScalarType> kinEngType;
    typedef mcpack::hamiltonian::leapfrog<potEngType,kinEngType> integratorType;
    typedef typename mcpack::utils::randomSTD<realScalarType> randVarGenType;
    typedef mcpack::hamiltonian::HMCProposal<integratorType,randVarGenType> HMCProposalType;

    const indexType N = 10;

    realVectorType mu = realVectorType::Zero(N);
    realMatrixType sigmaInv = realMatrixType::Identity(N,N);
    realMatrixType MInv = realMatrixType::Identity(N,N);
    realVectorType q0 = realVectorType::Random(N);

    const realScalarType eps = 1;
    const indexType Nsteps = 10;

    potEngType G(mu,sigmaInv);
    kinEngType K(MInv);

    integratorType Lp(G,K);

    HMCProposalType prop(Lp,eps,Nsteps);

    realScalarType dH = 0;
    randVarGenType RVGen(12346l);

    prop.propose(q0,dH,RVGen);
}
