/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_RANDOM_HPP
#define MCPACK_RANDOM_HPP

namespace mcpack{ namespace utils {

    template<typename _RealType>
    class RandomVariateGenerator
    {
    public:
        static_assert(std::is_floating_point<_RealType>::value,
            "Parameter should be a floating point type");
        
        typedef _RealType RealType;
        typedef _RealType InputType;
        typedef _RealType ResultType;
        typedef std::mt19937_64 EngineType;
        typedef EngineType::result_type SeedType;
        typedef std::normal_distribution<RealType> NormDistType;
        typedef std::uniform_real_distribution<RealType> UniDistType;
        
        RandomVariateGenerator(SeedType seed=0)
        {
            Seed(seed);
        }
        
        ~RandomVariateGenerator()
        {
        
        }
        
        void Seed(SeedType sd)
        {
            m_eng.seed(sd);
        }
        
        ResultType Normal()
        {
            return m_norm(m_eng);
        }

        ResultType Uniform()
        {
            return m_uni(m_eng);
        }

        void GetState(std::stringstream & state) const
        {
            state << m_eng;
        }

        void SetState(std::stringstream & state)
        {
            state >> m_eng ;
        }
        
    private:
        EngineType m_eng;
        NormDistType m_norm;
        UniDistType m_uni;
    };
    
}//namespace utils
}//namespace mcpack

#endif//MCPACK_BOOST_RANDOM_HPP
