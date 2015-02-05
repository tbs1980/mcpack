#ifndef MCPACK_HMCPROPOSAL_HPP
#define MCPACK_HMCPROPOSAL_HPP

namespace mcpack { namespace hamiltonian {

    template<class _discretisationType,class _randVarGenType>
    class hmcProposal
    {
    public:
        typedef _discretisationType discretisationType;
        typedef _randVarGenType randVarGenType;
        typedef typename discretisationType::realScalarType realScalarType;
        typedef typename discretisationType::indexType indexType;
        typedef typename discretisationType::realVectorType realVectorType;

        hmcProposal()
        :m_eps(0),m_NSteps(0)
        {

        }

        hmcProposal(discretisationType const & Discr,const realScalarType eps,const indexType NSteps)
        :m_Discr(Discr),m_eps(eps),m_NSteps(NSteps)
        {

        }

        void propose(realVectorType & q0,realScalarType & deltaH,randVarGenType & RVGen)
        {
            //need to check if the descretisation has been intialised
            const indexType NDim=q0.rows();
            realScalarType u=m_eps*RVGen.Uniform();
            realScalarType eps=m_eps*u;
            u=m_eps*RVGen.Uniform();
            indexType NSteps=(indexType)(u*(realScalarType)m_NSteps);

            realVectorType p0(NDim);
            for(indexType i=0;i<NDim;++i)
            {
                p0(i)=RVGen.Normal();
            }

            m_Discr.Integrate(q0,p0,eps,NSteps,deltaH);

        }


    private:
        discretisationType m_Discr;
        realScalarType m_eps;
        indexType m_NSteps;
    };

} //namespace hamiltonian
} //namespace mcpack


#endif //MCPACK_HMCPROPOSAL_HPP
