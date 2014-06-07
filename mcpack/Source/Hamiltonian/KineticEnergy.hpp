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


#ifndef MCPACK_KINETICENERGY_HPP
#define MCPACK_KINETICENERGY_HPP

namespace mcpack { namespace hamiltonian {

    template<typename _RealType>
    class GaussKineticEnergy
    {
    public:
        static_assert(std::is_floating_point<_RealType>::value,
            "PARAMETER SHOULD BE A FLOATING POINT TYPE");

        typedef _RealType RealType;
        typedef typename Eigen::Matrix<RealType, Eigen::Dynamic, 1> RealVectorType;
        typedef typename Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic> RealMatrixType;
        typedef typename RealVectorType::Index IndexType;
        typedef typename Eigen::LLT<RealMatrixType> LLTType;

        GaussKineticEnergy()
        :m_MInv(0,0)
        {

        }

        explicit GaussKineticEnergy(RealMatrixType const& MInv)
        :m_MInv(MInv),m_Chol(MInv.rows(),MInv.cols())
        {
            MCPACK_ASSERT(MInv.rows()==MInv.cols(),"Mass^-1 should be a square matrix: rows==cols");

            //find the inverse of the matrix
            LLTType lltOfMInv(m_MInv);

            MCPACK_ASSERT(lltOfMInv.info()==Eigen::Success,"Mass^-1 is not positive definite");

            m_Chol=lltOfMInv.matrixL();
        }
        
        void Evaluate(RealVectorType const & p, RealType & val,RealVectorType & dp) const
        {
            MCPACK_ASSERT(p.rows()==dp.rows(),"p and dp shoudl have the same dimensionality");
            dp=-m_MInv*p;
            val=0.5*p.transpose()*dp;
        }

        void Rotate(RealVectorType & p) const
        {
            p=m_Chol*p;
        }

        IndexType NDim(void) const
        {
            return m_MInv.rows();
        }

    private:
        RealMatrixType m_MInv;
        RealMatrixType m_Chol;
    };

    template<typename _RealType>
    class GaussKineticEnergyDiag
    {
    public:
        static_assert(std::is_floating_point<_RealType>::value,
            "PARAMETER SHOULD BE A FLOATING POINT TYPE");

        typedef _RealType RealType;
        typedef typename Eigen::Matrix<RealType, Eigen::Dynamic, 1> RealVectorType;
        typedef typename Eigen::Matrix<RealType, Eigen::Dynamic, 1> RealDiagMatrixType;
        typedef typename RealDiagMatrixType::Index IndexType;

        GaussKineticEnergyDiag()
        :m_MInv((IndexType)0)
        {

        }
        
        explicit GaussKineticEnergyDiag(RealDiagMatrixType const& MInv)
        :m_MInv(MInv)
        {}
        
        void Evaluate(RealVectorType const & p, RealType & val,RealVectorType & dp) const
        {
            MCPACK_ASSERT(p.rows()==dp.rows(),"p and dp shoudl have the same dimensionality");
            dp=-m_MInv.cwiseProduct(p);
            val=0.5*p.transpose()*dp;       
        }

        void Rotate(RealVectorType & p) const
        {
            for(IndexType i=0;i<p.rows();++i)
            {
                p(i)*=sqrt(m_MInv(i));
            }
        }

        IndexType NDim(void) const
        {
            return m_MInv.rows();
        }

    private:
        RealDiagMatrixType m_MInv;
    };


} //namespace hamiltonian
} //namespace mcpack

#endif //MCPACK_KINETICENERGY_HPP
