/* 
* 
* Copyright (C) 2014 Sreekumar Thaithara Balan
* 
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or (at
* your option) any later version.
* 
* This program is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

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
