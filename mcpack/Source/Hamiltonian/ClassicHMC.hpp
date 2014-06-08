/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_CLASSIC_HMC_HPP
#define MCPACK_CLASSIC_HMC_HPP

namespace mcpack { namespace hamiltonian {

    template<class _DiscretisationType>
    class ClassicHMC
    {
    public:
        typedef _DiscretisationType DiscretisationType;
        typedef typename DiscretisationType::RealType RealType;
        typedef typename DiscretisationType::IndexType IndexType;
        typedef typename DiscretisationType::RealVectorType RealVectorType;
        typedef typename Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic> RealMatrixType;
        typedef typename mcpack::utils::RandomVariateGenerator<RealType> RandVarGenType;
        typedef typename RandVarGenType::SeedType SeedType;

        ClassicHMC()
        :m_eps(0),m_NSteps(0),m_q0(0,0)
        {}

        ClassicHMC(DiscretisationType const & Discr,RealType eps,IndexType NSteps,
            SeedType seed,RealVectorType const& q0)
        :m_Discr(Discr),m_eps(eps),m_NSteps(NSteps),m_RVGen(seed),m_q0(q0)
        {
        }

        void Generate(RealMatrixType & Samples)
        {
            IndexType NSamples=Samples.rows();
            IndexType NDim=Samples.cols();

            MCPACK_ASSERT(m_Discr.NDim()==NDim,
                "DiscretisationType and Samples should have the same dimensionality.");

            IndexType iter=0;
            IndexType samp=0;

            while(samp < NSamples)
            {
                RealType u=m_eps*m_RVGen.Uniform();
                RealType eps=m_eps*u;
                u=m_eps*m_RVGen.Uniform();
                IndexType NSteps=(IndexType)(u*(RealType)m_NSteps);
                
                RealVectorType p0(NDim);
                for(IndexType i=0;i<NDim;++i)
                {
                    p0(i)=m_RVGen.Normal();
                }

                RealType dH=0;
                RealVectorType q1(m_q0);
                m_Discr.Integrate(q1,p0,eps,NSteps,dH);

                u=m_eps*m_RVGen.Uniform();
                if(u < exp(-dH))
                {
                    m_q0=q1;
                    Samples.row(samp)=m_q0;
                    ++samp;
                }
                
                iter++;
            }

            m_AccRate=(RealType)samp/(RealType)iter;
        }

        RealType GetAcceptanceRate(void) const
        {
            return m_AccRate;
        }

        void SetSeed(unsigned long seed)
        {
            m_RVGen.Seed(seed);
        }

        void GetRandState(std::stringstream & RNGstate) const
        {
            m_RVGen.GetState(RNGstate);
        }

        void SetRandState(std::stringstream & RNGstate)
        {
            m_RVGen.SetState(RNGstate);
        }

        void SetStartPoint(RealVectorType const & q0)
        {
            m_q0=q0;
        }

    private:
        DiscretisationType m_Discr;
        RealType m_eps;
        IndexType m_NSteps;
        RandVarGenType m_RVGen;
        RealVectorType m_q0;
        RealType m_AccRate;
    };


} //namespace hamiltonian
} //namespace mcpack


#endif //MCPACK_CLASSIC_HMC_HPP
