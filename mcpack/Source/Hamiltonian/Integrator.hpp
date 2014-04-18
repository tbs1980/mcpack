#ifndef MCPACK_INTEGRATOR_HPP
#define MCPACK_INTEGRATOR_HPP

namespace mcpack { namespace hamiltonian {

	template<class _PotEngType,class _KinEngType>
	class LeapFrog
	{
	public:
		typedef _PotEngType PotEngType;
		typedef _KinEngType KinEngType;

		static_assert(typeid(PotEngType::RealType)==typeid(KinEngType::RealType),
			"Parameter should be a floating point type");

		typedef typename PotEngType::RealType RealType;
		typedef typename Eigen::Matrix<RealType, Eigen::Dynamic, 1> RealVectorType;
		typedef typename RealVectorType::Index IndexType;

		LeapFrog(PotEngType const & pe,KinEngType const& ke,RealType eps,IndexType N)
		:m_pe(pe),m_ke(ke),m_eps(eps),m_N(N)
		{
			MCPACK_ASSERT(m_pe.NDim()==m_ke.NDim(),
				"Pot-Energy and Kin-Energy objects should have the same dimensionality.");
			MCPACK_ASSERT(m_eps>0 and m_eps <2,"For stability of the LeapFrog, we require 0<eps<2");
			MCPACK_ASSERT(m_N>0,"Number of LeapFrog steps should be a positive integer.");
		}
		
		~LeapFrog(){}

		LeapFrog(LeapFrog const & other)
		{
			m_pe=other.m_pe;
			m_ke=other.m_ke;
		}
		
		void Integrate(RealVectorType & x,RealVectorType & p) const
		{
			MCPACK_ASSERT(x.rows()==p.rows(),
				"position and momentum should have the same number of dimensions");
			
			const IndexType N=x.rows();
			RealVectorType q=RealVectorType::Zero(N);
			RealVectorType g=RealVectorType::Zero(N);
			RealType valp=0;
			RealType valk=0;
			m_pe.Evaluate(x,valp,g);
			m_ke.Evaluate(p,valk,q);

			for(IndexType i=0;i<N;++i)
			{
				RealVectorType pe2=p+0.5*g;
				m_ke.Evaluate(pe2,valk,q);
				x=x-q;
				m_pe.Evaluate(x,valp,g);
				p=pe2+0.5*g;				
			}
		}

		void SetEps(RealType eps)
		{
			MCPACK_ASSERT(eps>0 and eps <2,"For stability of the LeapFrog, we require 0<eps<2");
			m_eps=eps;
		}

		void SetN(IndexType N)
		{
			MCPACK_ASSERT(N>0,"Number of LeapFrog steps should be a positive integer.");
			m_N=N;
		}
		
	private:
		 PotEngType m_pe;
		 KinEngType m_ke;
		 RealType m_eps;
		 RealType m_N;
		
	};



} //namespace hamiltonian
} //namespace mcpack


#endif //MCPACK_INTEGRATOR_HPP