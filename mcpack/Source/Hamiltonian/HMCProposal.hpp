#ifndef MCPACK_HMCProposal_HPP
#define MCPACK_HMCProposal_HPP

namespace mcpack { namespace hamiltonian {

    /**
     * \ingroup Hamiltonian
     *
     * \class HMCProposal
     *
     * \brief A class that implemets a proposal for Hamiltonian sampling
     *
     * \tparam _discretisationType the discretisation of the Hamilotnian integral
     * \tparam _randVarGenType random variate generator type
     */
    template<class _discretisationType,class _randVarGenType>
    class HMCProposal
    {
    public:

        /**
         * \typedef _discretisationType discretisationType;
         * \brief the discretisation of the phase space integral
         */
        typedef _discretisationType discretisationType;

        /**
         * \typedef _randVarGenType randVarGenType;
         * \brief the random variate generator type
         */
        typedef _randVarGenType randVarGenType;

        /**
         * \typedef typename discretisationType::realScalarType realScalarType;
         * \brief real floating point type
         */
        typedef typename discretisationType::realScalarType realScalarType;

        /**
         * \typedef typename discretisationType::indexType indexType;
         * \brief the integral type
         */
        typedef typename discretisationType::indexType indexType;

        /**
         * \typedef typename discretisationType::realVectorType realVectorType;
         * \brief real vector type
         */
        typedef typename discretisationType::realVectorType realVectorType;

        static_assert(std::is_floating_point<realScalarType>::value,
            "PARAMETER SHOULD BE A FLOATING POINT TYPE");

        /**
         * \brief the default constructor
         */
        HMCProposal()
        :m_eps(0),m_numSteps(0)
        {

        }

        /**
         * \brief a constructor that sets up the HMC proposal
         *
         * \param discr discretisation, e.g. leapfrog
         * \param eps espilon or maximum value of the step-size factor
         * \param numSteps maximum number of steps allowed
         */
        HMCProposal(discretisationType const & discr,const realScalarType eps,const indexType numSteps)
        :m_discr(discr),m_eps(eps),m_numSteps(numSteps)
        {
            BOOST_ASSERT_MSG(eps>0 and eps <2,"For stability of the leapfrog, we require 0<eps<2");
            BOOST_ASSERT_MSG(numSteps>0 ,"Number of steps should be a positive integer");
        }

        /**
         * \brief A function that proposes a new position
         *
         * \param q0     at input the current position, at output the new psition
         * \param deltaH difference in the Hamilronians
         * \param rvGen  random variate generator
         *
         * Here we assume that the discretisation in general requires the (q,p)
         * position in phase space, a step-size eps and the number steps to
         * take. Examples of this kind of discretisation are Euler's method,
         * Modified Euler's method and leapfrog method.
         * TODO a full descrption using maths
         */
        void propose(realVectorType & q0,realScalarType & deltaH,
            randVarGenType & rvGen) const
        {
            // TODO need to check if the descretisation has been intialised

            const indexType numDims=q0.rows();

            // randomise the step-size and the number of steps

            // randomise the step-size
            realScalarType u = rvGen.uniform(); //was m_eps*rvGen.uniform();
            const realScalarType eps = m_eps*u;

            // randomise the number of steps
            // we can do a reverse scaling here, e.g. u = 1. - u
            u = rvGen.uniform(); // was u = m_eps*rvGen.uniform();
            const indexType numSteps = (indexType)(u*(realScalarType)m_numSteps);
            // TODO can we use unifrom_int() here?

            // generate a random momentum vector
            realVectorType p0(numDims);
            for(indexType i=0;i<numDims;++i)
            {
                p0(i) = rvGen.normal();
            }

            // integrate the phase space using discretisation
            m_discr.integrate(q0,p0,eps,numSteps,deltaH);
        }

        /**
         * \brief Retrun the number of dimensions of the discretisation type
         *
         * \return the number of diemensions of the discretisation
         */
        inline indexType numDims(void) const
        {
            return m_discr.numDims();
        }

    private:
        discretisationType m_discr;/**< discretisation type */
        realScalarType m_eps;/**< maximum step-size factor */
        indexType m_numSteps;/**< maximum number of steps */
    };

} //namespace hamiltonian
} //namespace mcpack


#endif //MCPACK_HMCProposal_HPP
