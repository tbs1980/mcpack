/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_RANDOM_HPP
#define MCPACK_RANDOM_HPP

namespace mcpack{ namespace utils {

    /**
     * \ingroup utils
     *
     * \class randomSTD
     *
     * \brief A class for generating random variates using c++11 methods.
     *
     * \tparam _realScalarType The floating point type
     *
     * This class genrates random variates using c++11 methods. The random
     * number engine is mt19937.
     */
    template<typename _realScalarType>
    class randomSTD
    {
    public:
        static_assert(std::is_floating_point<_realScalarType>::value,
        "_realScalarType should be a floating point type");

        /**
         * \typedef _realScalarType realScalarType
         * \brief floating point type
         */
        typedef _realScalarType realScalarType;

        /**
         * \typedef std::mt19937 rngType
         * \brief Mersenne twister random number generator type
         */
        typedef std::mt19937 rngType;

        /**
         * \typedef typename std::normal_distribution<realScalarType> normDistType
         * \brief Normal distribution type
         */
        typedef typename std::normal_distribution<realScalarType> normDistType;

        /**
         * \typedef typename std::uniform_real_distribution<realScalarType> uniRealDistType
         * \brief uniform real distribution type
         */
        typedef typename std::uniform_real_distribution<realScalarType> uniRealDistType;

        /**
         * \typedef rngType::result_type seedType
         * \brief random number generator seed type
         */
        typedef rngType::result_type seedType;

        /**
         * \typedef typename Eigen::Matrix<realScalarType,Eigen::Dynamic,1> realVectorType
         * \brief real vector type
         */
        typedef typename Eigen::Matrix<realScalarType,Eigen::Dynamic,1> realVectorType;

        /**
         * \typedef typename realVectorType::Index indexType
         * \brief integral type
         */
        typedef typename realVectorType::Index indexType;

        /**
         * \brief The default constructor
         * \param seed seed for the random number generator
         */
        explicit randomSTD(seedType seed = rngType::default_seed)
        :m_rng(seed)
        {
        }

        inline void seed(seedType const seed)
        {
            m_rng.seed(seed);
        }

        /**
         * \brief A function that returns Normal(0,1) random variate
         * \return Normal(0,1) random variate
         */
        inline realScalarType normal(void)
        {
            return m_normDist(m_rng);
        }

        /**
         * \brief A function that returns Uniform(0,1) random variate
         * \return Uniform(0,1) random variate
         */
        inline realScalarType uniform()
        {
            return m_uniRealDist(m_rng);
        }

        /**
         * \brief A function that returns a vector of Normal(0,1) random variate
         * \param x vector to which Normal(0,1) variates are to be returned
         */
        void normal(realVectorType & x)
        {
            for(indexType i=0;i<x.rows();++i)
            {
                x(i) = m_normDist(m_rng);
            }
        }

        /**
         * \brief get the random number state
         * \param state State of the random number generator to be returned
         */
        inline void getState(std::stringstream & state) const
        {
            state << m_rng;
        }

        /**
         * \brief set the random number state
         * \param state State of the random number generator to be set
         */
        inline void setState(std::stringstream  & state)// cannot use const & here
        {
            state >> m_rng ;
        }

    private:
        rngType m_rng; /**< random number generator engine */
        normDistType m_normDist; /**< normal distribution */
        uniRealDistType m_uniRealDist; /**< uniform real distribution */
    };

}//namespace utils
}//namespace mcpack

#endif//MCPACK_BOOST_RANDOM_HPP
