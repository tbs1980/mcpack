#ifndef MCPACK_HMCPROPOSAL_HPP
#define MCPACK_HMCPROPOSAL_HPP

namespace mcpack { namespace hamiltonian {

    template<class _DiscretisationType,class _RandVarGenType>
    class HMCProposal
    {
    public:
        typedef _DiscretisationType DiscretisationType;
        typedef _RandVarGenType RandVarGenType;
        typedef typename DiscretisationType::RealType RealType;
        typedef typename DiscretisationType::IndexType IndexType;
        typedef typename DiscretisationType::RealVectorType RealVectorType;

        HMCProposal()
        :m_eps(0),m_NSteps(0)
        {

        }

        HMCProposal(DiscretisationType const & Discr,const RealType eps,const IndexType NSteps)
        :m_Discr(Discr),m_eps(eps),m_NSteps(NSteps)
        {

        }

        void Propose(RealVectorType & q0,RealType & deltaH,RandVarGenType & RVGen)
        {
            //need to check if the descretisation has been intialised
            const IndexType NDim=q0.rows();
            RealType u=m_eps*RVGen.Uniform();
            RealType eps=m_eps*u;
            u=m_eps*RVGen.Uniform();
            IndexType NSteps=(IndexType)(u*(RealType)m_NSteps);

            RealVectorType p0(NDim);
            for(IndexType i=0;i<NDim;++i)
            {
                p0(i)=RVGen.Normal();
            }

            m_Discr.Integrate(q0,p0,eps,NSteps,deltaH);

        }


    private:
        DiscretisationType m_Discr;
        RealType m_eps;
        IndexType m_NSteps;
    };

} //namespace hamiltonian
} //namespace mcpack


#endif //MCPACK_HMCPROPOSAL_HPP