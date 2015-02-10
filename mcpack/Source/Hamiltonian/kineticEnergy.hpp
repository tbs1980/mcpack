/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#ifndef MCPACK_KINETICENERGY_HPP
#define MCPACK_KINETICENERGY_HPP

namespace mcpack { namespace hamiltonian {

    /**
     * \ingroup Hamiltonian
     *
     * \class gaussKineticEnergy
     *
     * \brief A class for computing Gaussain Kinetic Enegry.
     *
     * \tparam _realScalarType real floating point type
     */
    template<typename _realScalarType>
    class gaussKineticEnergy
    {
    public:
        static_assert(std::is_floating_point<_realScalarType>::value,
            "PARAMETER SHOULD BE A FLOATING POINT TYPE");

        /**
         * \typedef _realScalarType realScalarType;
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
         * \typedef typename realVectorType::Index indexType;
         * \brief integral type
         */
        typedef typename realVectorType::Index indexType;

        /**
         * \typedef typename Eigen::LLT<realMatrixType> LLTType;
         * \brief Cholesky decompostion type
         */
        typedef typename Eigen::LLT<realMatrixType> LLTType;

        /**
         * \brief The default constructor
         */
        gaussKineticEnergy()
        :m_mInv(0,0)
        {}

        /**
         * \brief The default constructor that allocates the memory.
         *
         * \param mInv the inverse of the kinetic energy matrix
         */
        explicit gaussKineticEnergy(realMatrixType const& mInv)
        :m_mInv(mInv),m_Chol(mInv.rows(),mInv.cols())
        {
            BOOST_ASSERT_MSG(mInv.rows()==mInv.cols(),"Mass^-1 should be a square matrix: rows==cols");

            //find the inverse of the mInv ie we need the mass matrix M
            LLTType lltOfmInv(m_mInv.inverse());

            BOOST_ASSERT_MSG(lltOfmInv.info()==Eigen::Success,"Mass matrix is not positive definite");

            m_Chol=lltOfmInv.matrixL();
        }

        /**
         * \brief evaluate the kinetic energy
         *
         * \param p momentum at which kinetic energy is to be calculated
         * \param val value kinetic energy at \a p
         * \param dp derivative of the kinetic energy at \a p
         */
        void evaluate(realVectorType const & p, realScalarType & val,realVectorType & dp) const
        {
            BOOST_ASSERT_MSG(p.rows()==dp.rows(),"p and dp shoudl have the same dimensionality");
            dp = -m_mInv*p;
            val = 0.5*p.transpose()*dp;
        }

        /**
         * \brief rotate the momentum \a p using the kinetic energy matrix
         *
         * \param p momentum at which kinetic energy is to be calculated
         */
        inline void rotate(realVectorType & p) const
        {
            p = m_Chol*p;
        }

        /**
         * \brief return the number of dimensions
         *
         * \return the number of dimensions of the kinetic energy matrix
         */
        inline indexType numDims(void) const
        {
            return m_mInv.rows();
        }

    private:
        realMatrixType m_mInv; /**< kinetic energy matrix */
        realMatrixType m_Chol; /**< Colesky decomposition of the kinetic energy matrix */
    };


    /**
    * \ingroup Hamiltonian
    *
    * \class gaussKineticEnergyDiag
    *
    * \brief A class for computing Gaussain Kinetic Enegry with diagonal energy matrix
    *
    * \tparam _realScalarType real floating point type
    */
    template<typename _realScalarType>
    class gaussKineticEnergyDiag
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
        gaussKineticEnergyDiag()
        :m_mInv((indexType)0)
        {

        }

        /**
        * \brief The default constructor that allocates the memory.
        *
        * \param mInv the inverse of the kinetic energy matrix
        */
        explicit gaussKineticEnergyDiag(realDiagMatrixType const& mInv)
        :m_mInv(mInv)
        {}

        /**
         * \brief evaluate the kinetic energy
         *
         * \param p momentum at which kinetic energy is to be calculated
         * \param val value kinetic energy at \a p
         * \param dp derivative of the kinetic energy at \a p
         */
        void evaluate(realVectorType const & p, realScalarType & val,realVectorType & dp) const
        {
            BOOST_ASSERT_MSG(p.rows()==dp.rows(),"p and dp shoudl have the same dimensionality");
            dp = -m_mInv.cwiseProduct(p);
            val = 0.5*p.transpose()*dp;
        }

        /**
         * \brief rotate the momentum \a p using the kinetic energy matrix
         *
         * \param p momentum at which kinetic energy is to be calculated
         */
        void rotate(realVectorType & p) const
        {
            for(indexType i=0;i<p.rows();++i)
            {
                p(i) *= sqrt(1./m_mInv(i));
            }
        }

        /**
         * \brief return the number of dimensions
         *
         * \return the number of dimensions of the kinetic energy matrix
         */
        indexType numDims(void) const
        {
            return m_mInv.rows();
        }

    private:
        realDiagMatrixType m_mInv; /**< diagonal kinetic energy matrix */
    };


} //namespace hamiltonian
} //namespace mcpack

#endif //MCPACK_KINETICENERGY_HPP
