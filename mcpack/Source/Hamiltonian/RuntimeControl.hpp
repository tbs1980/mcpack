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

		explicit RunCtrl_FiniteSamples(IndexType const MaxSamples,IndexType const PacketSize)
		:m_NumSamples(0),m_MaxSamples(MaxSamples),m_PacketSize(PacketSize)
		{
			MCPACK_ASSERT(MaxSamples>0,"Maximum number of samples should be a positive integer");
		}

		bool Continue() const
		{
			return m_NumSamples >= m_MaxSamples ? false : true;
		}

		void Save(MatrixType const & Samples)
		{
			IndexType n=Samples.rows();
			m_NumSamples+=n;
		}

	private:
		IndexType m_NumSamples;
		IndexType m_MaxSamples;
		IndexType m_PacketSize;
	};

}
}

#endif //MCPACK_RUNTIMECONTROL_HPP