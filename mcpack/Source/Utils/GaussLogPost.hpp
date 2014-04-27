#ifndef MCPACK_GAUSSLOGPOST_HPP
#define MCPACK_GAUSSLOGPOST_HPP

namespace mcpack { namespace utils {

	template<typename _RealType>
	class GaussPotentialEnergy
	{
	public:
		static_assert(std::is_floating_point<_RealType>::value,
			"PARAMETER SHOULD BE A FLOATING POINT TYPE");

		typedef _RealType RealType;
		typedef typename Eigen::Matrix<RealType, Eigen::Dynamic, 1> RealVectorType;
		typedef typename Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic> RealMatrixType;
		typedef typename RealMatrixType::Index IndexType;

		GaussPotentialEnergy()
		:m_mu(0,0),m_SigmaInv(0,0)
		{}

		GaussPotentialEnergy(RealMatrixType const& mu,RealMatrixType const& SigmaInv)
		:m_mu(mu),m_SigmaInv(SigmaInv)
		{
			MCPACK_ASSERT(SigmaInv.rows()==SigmaInv.cols(),
					"Sigma^-1 should be a square matrix: rows==cols");
			MCPACK_ASSERT(SigmaInv.rows()==mu.rows(),
					"Sigma^-1 and mu should have the same dimensionality");	
			MCPACK_ASSERT(m_mu.rows()>0,"mu should have at least one element");
		}

		void Evaluate(RealVectorType const & q, RealType & val,RealVectorType & dq) const
		{
			MCPACK_ASSERT(q.rows()==dq.rows() && q.cols()==dq.cols(),
					"q and dq shoudl have the same dimensionality");
			MCPACK_ASSERT(q.rows()==NDim(),"q and dq shoudl have the same dimensionality");

			dq=m_SigmaInv*(m_mu-q);
			val=-0.5*(m_mu-q).transpose()*dq;
		}

		IndexType NDim(void) const
		{
			return m_SigmaInv.rows();
		}

private:
	RealVectorType m_mu;
	RealMatrixType m_SigmaInv;
};


}//namespace utils
}//namespace mcpack

#endif //MCPACK_GAUSSLOGPOST_HPP
