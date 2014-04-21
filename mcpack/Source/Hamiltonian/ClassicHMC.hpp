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

#ifndef MCPACK_CLASSIC_HMC_HPP
#define MCPACK_CLASSIC_HMC_HPP

namespace mcpack { namespace hamiltonian {

	template<class _DiscretisationType>
	class ClassicHMC
	{
	public:
		typedef _DiscretisationType DiscretisationType;
		typedef typename DiscretisationType::RealType RealType;
		typedef typename DiscretisationType::IndexType IndexType;
		typedef typename DiscretisationType::RealVectorType RealVectorType;
		typedef typename Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic> RealMatrixType;
		typedef typename mcpack::utils::RandomVariateGenerator<RealType> RandVarGenType;
		typedef typename RandVarGenType::SeedType SeedType;

		ClassicHMC(DiscretisationType const & Discr,RealType eps,IndexType NSteps,SeedType seed)
		:m_Discr(Discr),m_eps(eps),m_NSteps(NSteps),m_RVGen(seed)
		{
		}

		void Generate(RealVectorType & q0,RealMatrixType & Samples,RealType & AccRate)
		{
			IndexType NSamples=Samples.rows();
			IndexType NDim=Samples.cols();

			MCPACK_ASSERT(m_Discr.NDim()==NDim,
				"DiscretisationType and Samples should have the same dimensionality.");

			IndexType iter=0;
			IndexType samp=0;

			while(samp < NSamples)
			{
				RealType u=m_eps*m_RVGen.Uniform();
				RealType eps=m_eps*u;
				u=m_eps*m_RVGen.Uniform();
				IndexType NSteps=(IndexType)(u*(RealType)m_NSteps);
				
				RealVectorType p0(NDim);
				for(IndexType i=0;i<NDim;++i)
				{
					p0(i)=m_RVGen.Normal();
				}

				RealType dH=0;
				RealVectorType q1(q0);
				m_Discr.Integrate(q1,p0,eps,NSteps,dH);

				u=m_eps*m_RVGen.Uniform();
				if(u < exp(-dH))
				{
					q0=q1;
					Samples.row(samp)=q0;
					++samp;
				}
				
				iter++;
			}

			AccRate=(RealType)samp/(RealType)iter;
		}

		void GetRandState(std::stringstream & RNGstate)
		{
			m_RVGen.GetState(RNGstate);
		}

	private:
		DiscretisationType m_Discr;
		RealType m_eps;
		IndexType m_NSteps;
		RandVarGenType m_RVGen;
	};


} //namespace hamiltonian
} //namespace mcpack


#endif //MCPACK_CLASSIC_HMC_HPP
