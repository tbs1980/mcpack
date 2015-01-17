/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_leapfrog_HPP
#define MCPACK_leapfrog_HPP

namespace mcpack { namespace hamiltonian {

    /**
     * \ingroup Hamiltonian
     *
     * \class leapfrog
     *
     * \brief A class that implemets leapfrog integrator
     *
     * \tparam _potEngType Potentail Energy type
     * \tparam _kinEngType Kinetic Energy type
     *
     * This class implements the leapfrog integrator. MORE INFO TO COME.
     */
    template<class _potEngType,class _kinEngType>
    class leapfrog
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

        /**
         *  \brief The default constructor.
         */
        leapfrog(){}

        /**
         * \brief A constructor that potential and kintetic enrgies as arguments
         *
         * \param G Potential energy
         * \param K Kinetic energy
         */
        leapfrog(potEngType const & G,kinEngType const& K)
        :m_G(G),m_K(K)
        {
            BOOST_ASSERT_MSG(m_G.numDims()==m_K.numDims(),
                "Pot-Energy and Kin-Energy objects should have the same dimensionality.");
        }

        /**
         * \brief Integrate the Hamiltonian
         *
         * \param q positon vector
         * \param p momentum vector
         * \param eps epsilon or the step size
         * \param numSteps number steps in the integration
         * \param deltaH deifference in Hamiltonian after integration
         *
         * This method integrates the Hamilotian from (p,q) to (p',q') through
         * numSteps steps and epsilon step size.
         */
        void integrate(realVectorType & q,realVectorType & p,const realType eps,
            const indexType numSteps,realType & deltaH) const
        {
            BOOST_ASSERT_MSG(q.rows()==p.rows(),
                "position and momentum should have the same number of dimensions");
            BOOST_ASSERT_MSG(q.rows()==m_G.numDims(),
                "position and momentum should have the same number of dimensions");

            BOOST_ASSERT_MSG(eps>0 and eps <2,"For stability of the leapfrog, we require 0<eps<2");

            if(numSteps==0) return;

            const indexType N=q.rows();

            m_K.rotate(p);

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

            m_K.evaluate(p,valK,dp);
            realType h1=-(valG+valK);
            deltaH=h1-h0;
        }

        /**
         * \brief Retrun the number of dimensions of the Hamiltonian
         *
         * \return the number of diemensions of the Hamiltonian
         */
        indexType numDims(void) const
        {
            return m_G.numDims();
        }

    private:
        potEngType m_G; /**< potential energy */
        kinEngType m_K; /**< kinetic energy */
    };

} //namespace hamiltonian
} //namespace mcpack


#endif //MCPACK_INTEGRATOR_HPP
