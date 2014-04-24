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

		RunCtrl_FiniteSamples(IndexType const NumSamples, IndexType const PacketSize,
			IndexType const NumBurn)
		:m_Samples(0),m_NumSamples(NumSamples),m_PacketSize(PacketSize),
			m_Burn(0),m_NumBurn(NumBurn)
		{
			MCPACK_ASSERT(m_NumSamples>0,"Maximum number of samples should be a positive integer");
			MCPACK_ASSERT(m_NumBurn>0,"Number of samples to be burned should be a positive integer");
			MCPACK_ASSERT(m_NumBurn<m_NumSamples,"NumBurn should be < MaxSamples");
		}

		bool Continue() const
		{
			return m_Samples >= m_NumSamples ? false : true;
		}

		void Save(MatrixType const & Samples)
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

	private:
		IndexType m_Samples;
		IndexType m_NumSamples;
		IndexType m_PacketSize;
		IndexType m_Burn;
		IndexType m_NumBurn;
	};

}
}

#endif //MCPACK_RUNTIMECONTROL_HPP