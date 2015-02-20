/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#define BOOST_TEST_MODULE sampler
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>

template<typename realScalarType>
void test_sampler_finite_samples(void)
{
    typedef mcpack::utils::gaussPotentialEnergy<realScalarType> potEngType;
    typedef typename potEngType::realVectorType realVectorType;
    typedef typename potEngType::realMatrixType realMatrixType;
    typedef typename realMatrixType::Index indexType;
    typedef mcpack::hamiltonian::gaussKineticEnergy<realScalarType> kinEngType;
    typedef mcpack::hamiltonian::leapfrog<potEngType,kinEngType> integratorType;
    typedef mcpack::utils::randomSTD<realScalarType> randVarGenType;
    typedef mcpack::hamiltonian::HMCProposal<integratorType,randVarGenType> HMCProposalType;
    typedef mcpack::hamiltonian::classicHMC<HMCProposalType> HMCType;
    typedef mcpack::hamiltonian::IOWriteAllParams<realMatrixType> IOType;
    typedef mcpack::hamiltonian::runCtrlFiniteSamples<realMatrixType> RCType;
    typedef mcpack::hamiltonian::sampler<HMCType,IOType,RCType> samplerType;


    //define the HMC
    const indexType N = 10;
    realVectorType mu = realVectorType::Zero(N);
    realMatrixType sigmaInv = realMatrixType::Identity(N,N);
    realMatrixType mInv = realMatrixType::Identity(N,N);
    realVectorType q0 = realVectorType::Random(N);

    const realScalarType eps=1;
    const indexType numsteps=10;

    potEngType G(mu,sigmaInv);
    kinEngType K(mInv);

    integratorType Lp(G,K);

    HMCProposalType prop(Lp,eps,numsteps);

    HMCType hmc(prop,q0,12346l);

    //define the IO
    const std::string fileName("testSampler.extract");
    IOType iowall(fileName);

    //define the Runtime Control
    const indexType numParams = 10;
    const indexType numSamples = 1000;
    const indexType numBurn = 2000;
    const indexType packetSize = 100;
    const std::string fileRoot("./testSampler");
    const bool silent = false;

    RCType runctrl(numParams,numSamples,packetSize,numBurn,fileRoot,silent);

    samplerType Smp(hmc,iowall,runctrl);

    Smp.run();
}

template<typename realScalarType>
void test_sampler_finite_samples_diag_KE_and_PE(void)
{
    typedef mcpack::utils::gaussPotentialEnergyDiag<realScalarType> potEngType;
    typedef typename potEngType::realVectorType realVectorType;
    typedef typename potEngType::realDiagMatrixType realDiagMatrixType;
    typedef typename realDiagMatrixType::Index indexType;
    typedef mcpack::hamiltonian::gaussKineticEnergyDiag<realScalarType> kinEngType;
    typedef mcpack::hamiltonian::leapfrog<potEngType,kinEngType> integratorType;
    typedef mcpack::utils::randomSTD<realScalarType> randVarGenType;
    typedef mcpack::hamiltonian::HMCProposal<integratorType,randVarGenType> HMCProposalType;
    typedef mcpack::hamiltonian::classicHMC<HMCProposalType> HMCType;
    typedef typename HMCType::realMatrixType realMatrixType;
    typedef mcpack::hamiltonian::IOWriteAllParams<realMatrixType> IOType;
    typedef mcpack::hamiltonian::runCtrlFiniteSamples<realMatrixType> RCType;
    typedef mcpack::hamiltonian::sampler<HMCType,IOType,RCType> samplerType;


    //define the HMC
    const indexType N = 10;
    realVectorType mu = realVectorType::Zero(N);
    realDiagMatrixType sigmaInv(N);
    realDiagMatrixType mInv(N);
    realVectorType q0 = realVectorType::Random(N);
    for(indexType i=0;i<N;++i)
    {
        sigmaInv(i) = 1;
        mInv(i) = 1;
    }

    const realScalarType eps = 1;
    const indexType numsteps = 10;

    potEngType G(mu,sigmaInv);
    kinEngType K(mInv);

    integratorType Lp(G,K);

    HMCProposalType prop(Lp,eps,numsteps);

    HMCType hmc(prop,q0,12346l);

    //define the IO
    const std::string fileName("testSamplerDiag.extract");
    IOType iowall(fileName);

    //define the Runtime Control
    const indexType numParams = 10;
    const indexType numSamples = 1000;
    const indexType numBurn = 2000;
    const indexType packetSize = 100;
    const std::string fileRoot("./testSamplerDiag");
    const bool silent = false;

    RCType runctrl(numParams,numSamples,packetSize,numBurn,fileRoot,silent);

    samplerType Smp(hmc,iowall,runctrl);

    Smp.run();
}

BOOST_AUTO_TEST_CASE(sampler_finite_samples)
{
    test_sampler_finite_samples<float>();
    test_sampler_finite_samples<double>();
    test_sampler_finite_samples<long double>();
}

BOOST_AUTO_TEST_CASE(sampler_finite_samples_diag_KE_and_PE)
{
    test_sampler_finite_samples_diag_KE_and_PE<float>();
    test_sampler_finite_samples_diag_KE_and_PE<double>();
    test_sampler_finite_samples_diag_KE_and_PE<long double>();
}
