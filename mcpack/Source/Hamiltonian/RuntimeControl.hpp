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

#ifndef MCPACK_RUNTIMECONTROL_HPP
#define MCPACK_RUNTIMECONTROL_HPP


namespace mcpack { namespace hamiltonian {


	template<class _MatrixType>
	class RunCtrl_FiniteSamples
	{
	public:
		typedef _MatrixType MatrixType;
		typedef typename MatrixType::Index IndexType;

		RunCtrl_FiniteSamples()
		:m_NumParas(0),m_Samples(0),m_NumSamples(0),
		m_PacketSize(0),m_Burn(0),m_NumBurn(0)
		{}

		RunCtrl_FiniteSamples(IndexType const NumParas,IndexType const NumSamples, 
			IndexType const PacketSize,	IndexType const NumBurn,std::string const& root)
		:m_NumParas(NumParas),m_Samples(0),m_NumSamples(NumSamples),
		m_PacketSize(PacketSize),m_Burn(0),m_NumBurn(NumBurn),m_root(root)
		{
			MCPACK_ASSERT(m_NumSamples>0,"Maximum number of samples should be a positive integer");
			MCPACK_ASSERT(m_NumBurn>0,"Number of samples to be burned should be a positive integer");
			MCPACK_ASSERT(m_NumBurn<m_NumSamples,"NumBurn should be < MaxSamples");
		}

		bool Continue() const
		{
			return m_Samples >= m_NumSamples ? false : true;
		}

		void Add(MatrixType const & Samples)
		{
			IndexType n=Samples.rows();
			if(m_Burn >= m_NumBurn)
			{
				m_Samples+=n;
			}
			else
			{
				m_Burn+=n;
			}
		}

		IndexType NumParas(void) const
		{
			return m_NumParas;
		}

		IndexType PacketSize(void) const
		{
			return m_PacketSize;
		}

		std::string Root(void) const
		{
			return m_root;
		}

	private:
		IndexType m_NumParas;
		IndexType m_Samples;
		IndexType m_NumSamples;
		IndexType m_PacketSize;
		IndexType m_Burn;
		IndexType m_NumBurn;
		std::string m_root;
	};

}
}

#endif //MCPACK_RUNTIMECONTROL_HPP