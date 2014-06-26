/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_CLASSIC_HMC_HPP
#define MCPACK_CLASSIC_HMC_HPP

namespace mcpack { namespace hamiltonian {

    template<class _ProposalType>
    class ClassicHMC
    {
    public:
        typedef _ProposalType ProposalType;
        typedef typename ProposalType::RealType RealType;
        typedef typename ProposalType::IndexType IndexType;
        typedef typename ProposalType::RealVectorType RealVectorType;
        typedef typename ProposalType::RandVarGenType RandVarGenType;
        typedef typename Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic> RealMatrixType;
        //typedef typename mcpack::utils::RandomVariateGenerator<RealType> RandVarGenType;
        typedef typename RandVarGenType::SeedType SeedType;

        ClassicHMC()
        :m_q0(0,0),m_AccRate(0),m_RVGen(0)
        {}

        ClassicHMC(ProposalType const & Prop,RealVectorType const& q0,const SeedType seed)
        :m_Prop(Prop),m_q0(q0),m_AccRate(0),m_RVGen(seed)
        {
        }

        void Generate(RealMatrixType & Samples)
        {
            IndexType NSamples=Samples.rows();
            IndexType NDim=Samples.cols();

            //MCPACK_ASSERT(m_Discr.NDim()==NDim,
                //"DiscretisationType and Samples should have the same dimensionality.");

            IndexType iter=0;
            IndexType samp=0;

            while(samp < NSamples)
            {
                RealType dH=0;
                RealVectorType q1(m_q0);
                m_Prop.Propose(q1,dH,m_RVGen);

                RealType u=m_RVGen.Uniform();
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

        //the folling methods are required for 
        //setting seeds from outside
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
        ProposalType m_Prop;
        RealVectorType m_q0;
        RealType m_AccRate;
        RandVarGenType m_RVGen;
    };


} //namespace hamiltonian
} //namespace mcpack


#endif //MCPACK_CLASSIC_HMC_HPP
