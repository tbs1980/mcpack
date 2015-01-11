/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_LEAPFROG_HPP
#define MCPACK_LEAPFROG_HPP

namespace mcpack { namespace hamiltonian {

    template<class _potEngType,class _kinEngType>
    class leapFrog
    {
    public:
        typedef _potEngType potEngType;
        typedef _kinEngType kinEngType;

        typedef typename potEngType::realType realType;
        typedef typename potEngType::realVectorType realVectorType;
        typedef typename realVectorType::Index indexType;

    private:
        typedef typename kinEngType::realType realTypeKE;
        typedef typename kinEngType::realVectorType realVectorTypeKE;

        static_assert(std::is_floating_point<realType>::value,
            "PARAMETER SHOULD BE A FLOATING POINT TYPE");

        static_assert(std::is_same<realType,realTypeKE>::value,
            "POTENTIAL ENERGY AND KINTETIC ENERGY SHOULD SHOULD HAVE THE SAME FLOATING POINT TYPE");

        static_assert(std::is_same<realVectorType,realVectorTypeKE>::value,
            "POTENTIAL ENERGY AND KINTETIC ENERGY SHOULD SHOULD HAVE THE SAME FLOATING POINT TYPE");
    public:

        leapFrog(){}

        leapFrog(potEngType const & G,kinEngType const& K)
        :m_G(G),m_K(K)
        {
            //MCPACK_ASSERT(m_G.numDim()==m_K.numDim(),
            //    "Pot-Energy and Kin-Energy objects should have the same dimensionality.");
        }

        void integrate(realVectorType & q,realVectorType & p,const realType eps,
            const indexType numSteps,realType & deltaH) const
        {
            //MCPACK_ASSERT(q.rows()==p.rows(),
            //    "position and momentum should have the same number of dimensions");
            //MCPACK_ASSERT(q.rows()==m_G.numDim(),
            //    "position and momentum should have the same number of dimensions");

            //MCPACK_ASSERT(eps>0 and eps <2,"For stability of the LeapFrog, we require 0<eps<2");

            if(numSteps==0) return;

            const indexType N=q.rows();

            m_K.Rotate(p);

            realVectorType dp=realVectorType::Zero(N);
            realVectorType dq=realVectorType::Zero(N);
            realType valG=0;
            realType valK=0;
            m_G.evaluate(q,valG,dq);
            m_K.evaluate(p,valK,dp);

            realType h0=-(valG+valK);

            //take half a step
            p=p+0.5*eps*dq;

            for(indexType i=0;i<numSteps;++i)
            {
                //now full steps
                m_K.evaluate(p,valK,dp);
                q=q-eps*dp;
                m_G.evaluate(q,valG,dq);
                p=p+eps*dq;
            }

            //move the momentum back half a step
            p=p-0.5*eps*dq;

            m_K.Evaluate(p,valK,dp);
            realType h1=-(valG+valK);
            deltaH=h1-h0;
        }

        indexType numDim(void) const
        {
            return m_G.numDim();
        }

    private:
        potEngType m_G;
        kinEngType m_K;
    };

} //namespace hamiltonian
} //namespace mcpack


#endif //MCPACK_INTEGRATOR_HPP
