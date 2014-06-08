/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_INTEGRATOR_HPP
#define MCPACK_INTEGRATOR_HPP

namespace mcpack { namespace hamiltonian {

    template<class _PotEngType,class _KinEngType>
    class LeapFrog
    {
    public:
        typedef _PotEngType PotEngType;
        typedef _KinEngType KinEngType;

        typedef typename PotEngType::RealType RealType;
        typedef typename PotEngType::RealVectorType RealVectorType;
        typedef typename RealVectorType::Index IndexType;

    private:
        typedef typename KinEngType::RealType RealTypeKE;
        typedef typename KinEngType::RealVectorType RealVectorTypeKE;

        static_assert(std::is_floating_point<RealType>::value,
            "PARAMETER SHOULD BE A FLOATING POINT TYPE");

        static_assert(std::is_same<RealType,RealTypeKE>::value,
            "POTENTIAL ENERGY AND KINTETIC ENERGY SHOULD SHOULD HAVE THE SAME FLOATING POINT TYPE");

        static_assert(std::is_same<RealVectorType,RealVectorTypeKE>::value,
            "POTENTIAL ENERGY AND KINTETIC ENERGY SHOULD SHOULD HAVE THE SAME FLOATING POINT TYPE");
    public:

        LeapFrog(){}

        LeapFrog(PotEngType const & G,KinEngType const& K)
        :m_G(G),m_K(K)
        {
            MCPACK_ASSERT(m_G.NDim()==m_K.NDim(),
                "Pot-Energy and Kin-Energy objects should have the same dimensionality.");
        }
        
        void Integrate(RealVectorType & q,RealVectorType & p,const RealType eps,
            const IndexType NSteps,RealType & deltaH) const
        {
            MCPACK_ASSERT(q.rows()==p.rows(),
                "position and momentum should have the same number of dimensions");
            MCPACK_ASSERT(q.rows()==m_G.NDim(),
                "position and momentum should have the same number of dimensions");

            MCPACK_ASSERT(eps>0 and eps <2,"For stability of the LeapFrog, we require 0<eps<2");
            
            if(NSteps==0) return;
            
            const IndexType N=q.rows();

            m_K.Rotate(p);

            RealVectorType dp=RealVectorType::Zero(N);
            RealVectorType dq=RealVectorType::Zero(N);
            RealType valG=0;
            RealType valK=0;
            m_G.Evaluate(q,valG,dq);
            m_K.Evaluate(p,valK,dp);

            RealType h0=-(valG+valK);

            //take half a step 
            p=p+0.5*eps*dq;

            for(IndexType i=0;i<NSteps;++i)
            {
                //now full steps
                m_K.Evaluate(p,valK,dp);
                q=q-eps*dp;
                m_G.Evaluate(q,valG,dq);
                p=p+eps*dq;             
            }
            
            //move the momentum back half a step
            p=p-0.5*eps*dq;
            
            m_K.Evaluate(p,valK,dp);
            RealType h1=-(valG+valK);
            deltaH=h1-h0;
        }

        IndexType NDim(void) const
        {
            return m_G.NDim();
        }

    private:
        PotEngType m_G;
        KinEngType m_K;
    };

} //namespace hamiltonian
} //namespace mcpack


#endif //MCPACK_INTEGRATOR_HPP
