/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#define BOOST_TEST_MODULE classicHMC
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <mcpack/CoreHeaders.hpp>

template<typename realScalarType>
void test_classic_hmc_10(void)
{
    typedef typename mcpack::utils::gaussPotentialEnergy<realScalarType> potEngType;
    typedef typename potEngType::realVectorType realVectorType;
    typedef typename potEngType::realMatrixType realMatrixType;
    typedef typename realMatrixType::Index indexType;
    typedef mcpack::hamiltonian::gaussKineticEnergy<realScalarType> kinEngType;
    typedef mcpack::hamiltonian::leapfrog<potEngType,kinEngType> integratorType;
    typedef mcpack::utils::randomSTD<realScalarType> randVarGenType;
    typedef mcpack::hamiltonian::HMCProposal<integratorType,randVarGenType> HMCProposalType;
    typedef mcpack::hamiltonian::classicHMC<HMCProposalType> HMCType;

    const indexType N=10;

    realVectorType mu=realVectorType::Zero(N);
    realMatrixType sigmaInv=realMatrixType::Identity(N,N);
    realMatrixType mInv=realMatrixType::Identity(N,N);
    realVectorType q0=realVectorType::Random(N);

    const realScalarType eps=1;
    const indexType Nsteps=10;

    potEngType G(mu,sigmaInv);
    kinEngType K(mInv);

    integratorType Lp(G,K);

    HMCProposalType prop(Lp,eps,Nsteps);

    HMCType hmc(prop,q0,12346l);

    const indexType numsamples=1000;

    realMatrixType samples(numsamples,N);

    hmc.generate(samples);

    std::cout<<"Acceptace Rate= "<<hmc.getAcceptanceRate()<<std::endl;

    BOOST_REQUIRE(0.91 < hmc.getAcceptanceRate() and hmc.getAcceptanceRate() < 0.95);

}

template<typename realScalarType>
void test_classic_hmc_100(void)
{
    typedef typename mcpack::utils::gaussPotentialEnergy<realScalarType> potEngType;
    typedef typename potEngType::realVectorType realVectorType;
    typedef typename potEngType::realMatrixType realMatrixType;
    typedef typename realMatrixType::Index indexType;
    typedef mcpack::hamiltonian::gaussKineticEnergy<realScalarType> kinEngType;
    typedef mcpack::hamiltonian::leapfrog<potEngType,kinEngType> integratorType;
    typedef mcpack::utils::randomSTD<realScalarType> randVarGenType;
    typedef mcpack::hamiltonian::HMCProposal<integratorType,randVarGenType> HMCProposalType;
    typedef mcpack::hamiltonian::classicHMC<HMCProposalType> HMCType;

    const indexType N=100;

    realVectorType mu=realVectorType::Zero(N);
    realMatrixType sigmaInv=realMatrixType::Identity(N,N);
    realMatrixType mInv=realMatrixType::Identity(N,N);
    realVectorType q0=realVectorType::Random(N);

    const realScalarType eps=1;
    const indexType Nsteps=10;

    potEngType G(mu,sigmaInv);
    kinEngType K(mInv);

    integratorType Lp(G,K);

    HMCProposalType prop(Lp,eps,Nsteps);

    HMCType hmc(prop,q0,12346l);

    const indexType numsamples=1000;

    realMatrixType samples(numsamples,N);

    hmc.generate(samples);

    std::cout<<"Acceptace Rate= "<<hmc.getAcceptanceRate()<<std::endl;

    BOOST_REQUIRE(0.77 < hmc.getAcceptanceRate() and hmc.getAcceptanceRate() < 0.81);

}

template<typename realScalarType>
void test_classic_hmc_diag_mass_mat(void)
{
    typedef typename mcpack::utils::gaussPotentialEnergy<realScalarType> potEngType;
    typedef typename potEngType::realVectorType realVectorType;
    typedef typename potEngType::realMatrixType realMatrixType;
    typedef typename realMatrixType::Index indexType;
    typedef mcpack::hamiltonian::gaussKineticEnergyDiag<realScalarType> kinEngType;
    typedef typename kinEngType::realDiagMatrixType realDiagMatrixType;
    typedef mcpack::hamiltonian::leapfrog<potEngType,kinEngType> integratorType;
    typedef mcpack::utils::randomSTD<realScalarType> randVarGenType;
    typedef mcpack::hamiltonian::HMCProposal<integratorType,randVarGenType> HMCProposalType;
    typedef mcpack::hamiltonian::classicHMC<HMCProposalType> HMCType;

    const indexType N=100;

    realVectorType mu=realVectorType::Zero(N);
    realMatrixType sigmaInv=realMatrixType::Identity(N,N);
    realVectorType q0=realVectorType::Random(N);

    const realScalarType eps=1;
    const indexType Nsteps=10;

    realDiagMatrixType mInv(N);
    for(indexType i=0;i<N;++i)
    {
        mInv(i)=1;
    }

    potEngType G(mu,sigmaInv);
    kinEngType K(mInv);

    integratorType Lp(G,K);

    HMCProposalType prop(Lp,eps,Nsteps);

    HMCType hmc(prop,q0,12346l);

    const indexType numsamples=1000;

    realMatrixType samples(numsamples,N);

    hmc.generate(samples);

    std::cout<<"Acceptace Rate= "<<hmc.getAcceptanceRate()<<std::endl;

    BOOST_REQUIRE(0.77 < hmc.getAcceptanceRate() and hmc.getAcceptanceRate() < 0.82);

}

template<typename realScalarType>
void test_classic_hmc_diag_mass_mat_sigma_inv(void)
{
    typedef typename mcpack::utils::gaussPotentialEnergyDiag<realScalarType> potEngType;
    typedef typename potEngType::realVectorType realVectorType;
    typedef typename potEngType::realDiagMatrixType realDiagMatrixType;
    typedef typename realDiagMatrixType::Index indexType;
    typedef mcpack::hamiltonian::gaussKineticEnergyDiag<realScalarType> kinEngType;
    typedef mcpack::hamiltonian::leapfrog<potEngType,kinEngType> integratorType;
    typedef mcpack::utils::randomSTD<realScalarType> randVarGenType;
    typedef mcpack::hamiltonian::HMCProposal<integratorType,randVarGenType> HMCProposalType;
    typedef mcpack::hamiltonian::classicHMC<HMCProposalType> HMCType;
    typedef typename HMCType::realMatrixType realMatrixType;

    const indexType N=1000;

    realVectorType mu=realVectorType::Zero(N);
    realVectorType q0=realVectorType::Random(N);

    const realScalarType eps=1;
    const indexType Nsteps=10;

    realDiagMatrixType mInv(N);
    realDiagMatrixType sigmaInv(N);
   for(indexType i=0;i<N;++i)
    {
        realScalarType s=2;
        sigmaInv(i)=s;
        mInv(i)=1./s;
    }

    potEngType G(mu,sigmaInv);
    kinEngType K(mInv);

    integratorType Lp(G,K);

    HMCProposalType prop(Lp,eps,Nsteps);

    HMCType hmc(prop,q0,12346l);

    const indexType numsamples=1000;

    realMatrixType samples(numsamples,N);

    hmc.generate(samples);

    std::cout<<"Acceptace Rate= "<<hmc.getAcceptanceRate()<<std::endl;

    BOOST_REQUIRE(0.56 < hmc.getAcceptanceRate() and hmc.getAcceptanceRate() < 0.59);

}

template<typename realScalarType>
void test_classic_hmc_2d_correlated(void)
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

    const indexType N=2;

    realVectorType mu=realVectorType::Zero(N);
    realMatrixType sigmaInv=realMatrixType::Identity(N,N);
    realMatrixType mInv=realMatrixType::Identity(N,N);
    realVectorType q0=realVectorType::Random(N);


    //assign the off diagonal terms
    sigmaInv(0,1) = 0.5;
    sigmaInv(1,0) = 0.5;
    //set mInv = (sigmaInv)^-1
    mInv = sigmaInv.inverse();


    std::cout<<"sigmaInv \n"<<sigmaInv<<"\n"<<std::endl;
    std::cout<<"mInv \n"<<mInv<<"\n"<<std::endl;


    const realScalarType eps=1;
    const indexType Nsteps=10;

    potEngType G(mu,sigmaInv);
    kinEngType K(mInv);

    integratorType Lp(G,K);

    HMCProposalType prop(Lp,eps,Nsteps);

    HMCType hmc(prop,q0,12346l);

    const indexType numsamples=1000;

    realMatrixType samples(numsamples,N);

    hmc.generate(samples);

    std::cout<<"Acceptace Rate= "<<hmc.getAcceptanceRate()<<std::endl;

    BOOST_REQUIRE(0.96 < hmc.getAcceptanceRate() and hmc.getAcceptanceRate() < 0.98);
}

template<typename realScalarType>
void test_classic_hmc_10d_correlated(void)
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

    const indexType N=10;

    realVectorType mu=realVectorType::Zero(N);
    realVectorType q0=realVectorType::Random(N);

    realMatrixType Mat=realMatrixType::Random(N,N);
    realMatrixType Sigma=Mat*Mat.transpose();;
    realMatrixType sigmaInv=Sigma.inverse();
    realMatrixType mInv=Sigma;

    std::cout<<"sigmaInv \n"<<sigmaInv<<"\n"<<std::endl;
    std::cout<<"mInv \n"<<mInv<<"\n"<<std::endl;


    const realScalarType eps=1;
    const indexType Nsteps=10;

    potEngType G(mu,sigmaInv);
    kinEngType K(mInv);

    integratorType Lp(G,K);

    HMCProposalType prop(Lp,eps,Nsteps);

    HMCType hmc(prop,q0,12346l);

    const indexType numsamples=1000;

    realMatrixType samples(numsamples,N);

    hmc.generate(samples);

    std::cout<<"Acceptace Rate= "<<hmc.getAcceptanceRate()<<std::endl;

    BOOST_REQUIRE(0.93 < hmc.getAcceptanceRate() and hmc.getAcceptanceRate() < 0.94);
}

BOOST_AUTO_TEST_CASE(classic_hmc_10)
{
    test_classic_hmc_10<float>();
    test_classic_hmc_10<double>();
    test_classic_hmc_10<long double>();
}

BOOST_AUTO_TEST_CASE(classic_hmc_100)
{
    test_classic_hmc_100<float>();
    test_classic_hmc_100<double>();
    test_classic_hmc_100<long double>();
}

BOOST_AUTO_TEST_CASE(classic_hmc_diag_mass_mat)
{
    test_classic_hmc_diag_mass_mat<float>();
    test_classic_hmc_diag_mass_mat<double>();
    test_classic_hmc_diag_mass_mat<long double>();
}

BOOST_AUTO_TEST_CASE(classic_hmc_diag_mass_mat_sigma_inv)
{
    test_classic_hmc_diag_mass_mat_sigma_inv<float>();
    test_classic_hmc_diag_mass_mat_sigma_inv<double>();
    test_classic_hmc_diag_mass_mat_sigma_inv<long double>();
}

BOOST_AUTO_TEST_CASE(classic_hmc_2d_correlated)
{
    test_classic_hmc_2d_correlated<float>();
    test_classic_hmc_2d_correlated<double>();
    test_classic_hmc_2d_correlated<long double>();
}

BOOST_AUTO_TEST_CASE(classic_hmc_10d_correlated)
{
    test_classic_hmc_10d_correlated<float>();
    test_classic_hmc_10d_correlated<double>();
    test_classic_hmc_10d_correlated<long double>();
}
