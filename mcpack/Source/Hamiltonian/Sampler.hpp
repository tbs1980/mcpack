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

#ifndef MCPACK_SAMPLER_HPP
#define MCPACK_SAMPLER_HPP

namespace mcpack{ namespace hamiltonian{

	template<class _Engine,class _IOType,class _RunCtrlType>
	class Sampler
	{
	public:
		typedef _Engine EngineType;
		typedef _IOType IOType;
		typedef _RunCtrlType RunCtrlType;

		typedef typename EngineType::RealMatrixType RealMatrixType;
		
		Sampler(EngineType const & Eng,IOType const& IO,RunCtrlType const& RunCtrl)
		:m_Eng(Eng),m_IO(IO),m_RunCtrl(RunCtrl)
		{
			//m_RunCtrl.Load();
		}

		void Run()
		{
			RealMatrixType Samples(m_RunCtrl.PacketSize(),m_RunCtrl.NumParas());
			std::stringstream RandState;

			while(m_RunCtrl.Continue())
			{
				m_Eng.Generate(Samples);
				m_Eng.GetRandState(RandState);
				SaveRandState(RandState);
				m_RunCtrl.Add(Samples);
				m_IO.Write(Samples);
			}		
		}

		void SaveRandState(std::stringstream const & rs) const
		{
			std::string RngFileName=m_RunCtrl.Root()+std::string(".radnom");
			std::ofstream RngFile;
			RngFile.open(RngFileName.c_str(),std::ios::trunc);
			RngFile<<rs.str()<<std::endl;
			RngFile.close();
		}
		
	private:
		EngineType m_Eng;
		IOType m_IO;
		RunCtrlType m_RunCtrl;
	};

}//namespace hamiltonian
}//namespace mcpack

#endif //MCPACK_SAMPLER_HPP
