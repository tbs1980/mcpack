/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_GAUSSLOGPOST_HPP
#define MCPACK_GAUSSLOGPOST_HPP

namespace mcpack { namespace utils {

    /**
     * \ingroup Utils
     *
     * \class gaussPotentialEnergy
     *
     * \brief A class for computing Gaussian Potential Energy
     *
     * \tparam _realScalarType real floating point type
     */
    template<typename _realScalarType>
    class gaussPotentialEnergy
    {
    public:
        static_assert(std::is_floating_point<_realScalarType>::value,
            "PARAMETER SHOULD BE A FLOATING POINT TYPE");

        /**
         * \typedef _realScalarType realScalarType
         * \brief the floating point type
         */
        typedef _realScalarType realScalarType;

        /**
         * \typedef typename Eigen::Matrix<realScalarType, Eigen::Dynamic, 1> realVectorType;
         * \brief real vector type
         */
        typedef typename Eigen::Matrix<realScalarType, Eigen::Dynamic, 1> realVectorType;

        /**
         * \typedef typename Eigen::Matrix<realScalarType, Eigen::Dynamic, Eigen::Dynamic> realMatrixType;
         * \brief real matrix type
         */
        typedef typename Eigen::Matrix<realScalarType, Eigen::Dynamic, Eigen::Dynamic> realMatrixType;

        /**
         * \typedef typename realMatrixType::Index indexType;
         * \brief integral type
         */
        typedef typename realMatrixType::Index indexType;

        /**
         * \brief The default constructor
         */
        gaussPotentialEnergy()
        :m_mu(0,0),m_sigmaInv(0,0)
        {}

        /**
         * \brief A constructor that allocates the memory
         *
         * \param mu the mean of the Gaussian distribution
         * \param sigmaInv Inverse of the potential energy matrix
         */
        gaussPotentialEnergy(realMatrixType const& mu,realMatrixType const& sigmaInv)
        :m_mu(mu),m_sigmaInv(sigmaInv)
        {
            BOOST_ASSERT_MSG(sigmaInv.rows()==sigmaInv.cols(),
                    "Sigma^-1 should be a square matrix: rows==cols");
            BOOST_ASSERT_MSG(sigmaInv.rows()==mu.rows(),
                    "Sigma^-1 and mu should have the same dimensionality");
            BOOST_ASSERT_MSG(m_mu.rows()>0,"mu should have at least one element");
        }

        /**
         * \brief compute the potential energy
         *
         * \param q the position at which potential energy is to be calculated
         * \param val the potential energy at \a q
         * \param dq the derivative of the potential energy at \a q
         */
        void evaluate(realVectorType const & q, realScalarType & val,realVectorType & dq) const
        {
            BOOST_ASSERT_MSG(q.rows()==dq.rows() && q.cols()==dq.cols(),
                    "q and dq shoudl have the same dimensionality");
            BOOST_ASSERT_MSG(q.rows()==numDims(),"q and dq shoudl have the same dimensionality");

            dq = m_sigmaInv*(m_mu-q);
            val = -0.5*(m_mu-q).transpose()*dq;
        }

        /**
         * \brief return the number of dimensions of the potential energy matrix
         *
         * \return the number of dimensions of the potential energy matrix
         */
        inline indexType numDims(void) const
        {
            return m_sigmaInv.rows();
        }

    private:
        realVectorType m_mu; /**< the mean of the potential energy distribution */
        realMatrixType m_sigmaInv; /**< the inverse of the potential energy matrix */
    };


    /**
     * \ingroup Utils
     *
     * \class gaussPotentialEnergy
     *
     * \brief A class for computing Gaussian Potential Energy with diagonal matrix
     *
     * \tparam _realScalarType real floating point type
     */
    template<typename _realScalarType>
    class gaussPotentialEnergyDiag
    {
    public:
        static_assert(std::is_floating_point<_realScalarType>::value,
            "PARAMETER SHOULD BE A FLOATING POINT TYPE");

        typedef _realScalarType realScalarType;
        typedef typename Eigen::Matrix<realScalarType, Eigen::Dynamic, 1> realVectorType;
        typedef typename Eigen::Matrix<realScalarType, Eigen::Dynamic, 1> realDiagMatrixType;
        typedef typename realDiagMatrixType::Index indexType;

        /**
         * \brief The default constructor
         */
        gaussPotentialEnergyDiag()
        :m_mu(0,0),m_sigmaInv(0,0)
        {}

        /**
         * \brief A constructor that allocates the memory
         *
         * \param mu the mean of the Gaussian distribution
         * \param sigmaInv Inverse of the potential energy matrix
         */
        gaussPotentialEnergyDiag(realVectorType const& mu,realDiagMatrixType const& sigmaInv)
        :m_mu(mu),m_sigmaInv(sigmaInv)
        {
            BOOST_ASSERT_MSG(sigmaInv.rows()==mu.rows(),
                    "Sigma^-1 and mu should have the same dimensionality");
            BOOST_ASSERT_MSG(m_mu.rows()>0,"mu should have at least one element");
        }

        /**
         * \brief compute the potential energy
         *
         * \param q the position at which potential energy is to be calculated
         * \param val the potential energy at \a q
         * \param dq the derivative of the potential energy at \a q
         */
        void evaluate(realVectorType const & q, realScalarType & val,realVectorType & dq) const
        {
            BOOST_ASSERT_MSG(q.rows()==dq.rows() && q.cols()==dq.cols(),
                    "q and dq shoudl have the same dimensionality");
            BOOST_ASSERT_MSG(q.rows()==numDims(),"q and dq shoudl have the same dimensionality");

            dq = m_sigmaInv.cwiseProduct(m_mu-q);
            val = -0.5*(m_mu-q).transpose()*dq;
        }

        /**
         * \brief return the number of dimensions of the potential energy matrix
         *
         * \return the number of dimensions of the potential energy matrix
         */
        inline indexType numDims(void) const
        {
            return m_sigmaInv.rows();
        }

    private:

        realVectorType m_mu; /**< the mean of the potential energy distribution */
        realDiagMatrixType m_sigmaInv; /**< the inverse of the potential energy matrix */
    };


}//namespace utils
}//namespace mcpack

#endif //MCPACK_GAUSSLOGPOST_HPP
