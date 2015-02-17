/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_CLASSIC_HMC_HPP
#define MCPACK_CLASSIC_HMC_HPP

namespace mcpack { namespace hamiltonian {

    /**
     * \ingroup Hamiltonian
     *
     * \class HMCProposal
     *
     * \brief A class that implemets classic Hamiltonian sampling
     *
     * \tparam _proposalType Hamiltonian proposal type
     *
     * This class implements the classic Hamiltonian sampling using a
     * Hamiltonian proposal policy. The sampler is constructed with a
     * starting point and can return the requested number of samples.
     */
    template<class _proposalType>
    class classicHMC
    {
    public:

        /**
         * \typedef _proposalType proposalType;
         * \brief HMC proposal type
         */
        typedef _proposalType proposalType;

        /**
         * \typedef typename proposalType::realScalarType realScalarType;
         * \brief the floating point type
         */
        typedef typename proposalType::realScalarType realScalarType;

        /**
         * \typedef typename proposalType::indexType indexType;
         * \brief the integral type
         */
        typedef typename proposalType::indexType indexType;

        /**
         * \typedef typename proposalType::realVectorType realVectorType;
         * \brief real vector type
         */
        typedef typename proposalType::realVectorType realVectorType;

        /**
         * \typedef typename proposalType::randVarGenType randVarGenType;
         * \brief random variate generator type
         */
        typedef typename proposalType::randVarGenType randVarGenType;

        /**
         * \typedef typename Eigen::Matrix<realScalarType, Eigen::Dynamic, Eigen::Dynamic> realMatrixType;
         * \brief real matrix type
         */
        typedef typename Eigen::Matrix<realScalarType, Eigen::Dynamic, Eigen::Dynamic> realMatrixType;

        /**
         * \typedef typename randVarGenType::seedType seedType;
         * \brief random seed type
         */
        typedef typename randVarGenType::seedType seedType;

        static_assert(std::is_floating_point<realScalarType>::value,
            "PARAMETER SHOULD BE A FLOATING POINT TYPE");

        /**
         * \brief The default constructor
         */
        classicHMC()
        :m_q0(0,0),m_accRate(0)
        {
            // random variate generator is initiated with
            // its own random seed
        }

        /**
         * \brief A constructor that sets up the classic sampler
         *
         * \param prop a Hamiltonian proposal type
         * \param q0 the starting point
         * \param seed the random number generator seed
         */
        classicHMC(proposalType const & prop,realVectorType const& q0,const seedType seed)
        :m_prop(prop),m_q0(q0),m_accRate(0),m_RVGen(seed)
        {
        }

        /**
         * \brief Generate Monte Carlo samples
         * \param samples a container to hold the samples
         *
         * This method returns requested number of samples into the
         * real matrix \a samples. Each row in this matrix upon compeltion
         * will have a new samples. The matrix will have the dimession of
         * numSamples x numDims.
         * TODO we also need to output the log likelihood
         */
        void generate(realMatrixType & samples)
        {
            indexType numSamples = samples.rows();
            indexType numDims = samples.cols();

            BOOST_ASSERT_MSG(m_prop.numDims() == numDims,
                "DiscretisationType and samples should have the same dimensionality.");

            indexType iter = 0;
            indexType samp = 0;

            while(samp < numSamples)
            {
                realScalarType dH = 0;
                realVectorType q1(m_q0);
                m_prop.propose(q1,dH,m_RVGen);

                realScalarType u = m_RVGen.uniform();
                if(u < exp(-dH))
                {
                    m_q0=q1;
                    samples.row(samp)=m_q0;
                    ++samp;
                }
                iter++;
            }
            m_accRate = (realScalarType)samp/(realScalarType)iter;
        }

        /**
         * \brief return the acceptance reate
         * \return the acceptance rate of the samples
         */
        inline realScalarType getAcceptanceRate(void) const
        {
            return m_accRate;
        }

        /**
         * \brief set the seed of the random number generator
         * \param seed the random number seed to be set
         */
        inline void setSeed(seedType const seed)
        {
            m_RVGen.seed(seed);
        }

        /**
         * \brief extract the state of the random number generator
         * \param RNGstate the random number state to be returned
         */
        inline void getRandState(std::stringstream & RNGstate) const
        {
            m_RVGen.getState(RNGstate);
        }

        /**
         * \brief set the state of the random number generator
         * \param RNGstate random number generator state
         */
        inline void setRandState(std::stringstream & RNGstate)
        {
            m_RVGen.getState(RNGstate);
        }

        /**
         * \brief set a starting point to the classic samples
         * \param q0 the stating point
         */
        inline void setStartPoint(realVectorType const & q0)
        {
            m_q0=q0;
        }

    private:
        proposalType m_prop; /**< the Hamiltonian proposal */
        realVectorType m_q0; /**< starting point */
        realScalarType m_accRate; /**< acceptance rate */
        randVarGenType m_RVGen; /**< random variate generator */
    };


} //namespace hamiltonian
} //namespace mcpack


#endif //MCPACK_CLASSIC_HMC_HPP
