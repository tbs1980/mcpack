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
		typedef typename EngineType::RealVectorType RealVectorType;
		typedef typename EngineType::RealVectorType RealType;
		typedef typename EngineType::RealVectorType IndexType;
		
		Sampler(EngineType const & Eng,IOType const& IO,RunCtrlType const& RunCtrl)
		:m_Eng(Eng),m_IO(IO),m_RunCtrl(RunCtrl)
		{
			m_RngFileName=m_RunCtrl.Root()+std::string(".rng");
			m_LogFileName=m_RunCtrl.Root()+std::string(".log");
			m_BegFileName=m_RunCtrl.Root()+std::string(".beg");

			//first check if we have resume files
			std::ifstream RngFile;
			RngFile.open(m_RngFileName.c_str());
			std::ifstream LogFile;
			LogFile.open(m_LogFileName.c_str());
			std::ifstream BegFile;
			BegFile.open(m_BegFileName.c_str());

			if(RngFile.is_open() && LogFile.is_open() && BegFile.is_open())
			{
				std::cout<<"we have resume files\n"<<std::endl;
				if(m_RunCtrl.AreWeSampling(LogFile))
				{
					std::cout<<"we are sampling now\n"<<std::endl;
				}

				RngFile.close();
				BegFile.close();
				LogFile.close();
			}
			else
			{
				std::cout<<"We dont have resume files\n"<<std::endl;
			}


			//open the log file
			m_LogFile.open(m_LogFileName.c_str(),std::ios::app);
		}

		void Run()
		{
			RealMatrixType Samples(m_RunCtrl.PacketSize(),m_RunCtrl.NumParas());			
			
			while(m_RunCtrl.Continue())
			{
				//define the local stuff to logging
				std::stringstream RandState;
				std::stringstream LogState;
				RealVectorType ChainState;

				//generate samples
				m_Eng.Generate(Samples);
				m_RunCtrl.Add(Samples);

				//get states
				m_Eng.GetRandState(RandState);
				m_Eng.GetChainState(ChainState);

				//save states
				SaveRandState(RandState);
				SaveChainState(ChainState);

				//write log files
				m_RunCtrl.WriteLog(LogState);
				WriteLog(LogState);
				m_LogFile<<LogState.str()<<std::endl;

				//write console output
				if(!m_RunCtrl.Silent())
				{
					std::cout<<LogState.str()<<std::endl;
				}				

				//write extract
				m_IO.Write(Samples);
			}

			if(!m_RunCtrl.Silent())
			{
				std::cout<<"Sampling FINISHED\n"<<std::endl;				
			}

			m_LogFile<<"Sampling FINISHED\n"<<std::endl;
			m_LogFile.close();		
		}

	private:

		void SaveRandState(std::stringstream const & rs) const
		{
			std::ofstream RngFile;
			RngFile.open(m_RngFileName.c_str(),std::ios::trunc);
			RngFile<<rs.str()<<std::endl;
			RngFile.close();
		}

		void SaveChainState(RealVectorType const & cs) const
		{
			mcpack::utils::WriteMatrix2TextFile(cs,m_BegFileName,12,",");
		}

		void WriteLog(std::stringstream & LogHandle) const
		{
			LogHandle<<"Acceptance Rate= "<<m_Eng.GetAcceptanceRate()<<"\n\n"<<std::endl;
		}

		EngineType m_Eng;
		IOType m_IO;
		RunCtrlType m_RunCtrl;
		std::string m_RngFileName;
		std::string m_LogFileName;
		std::string m_BegFileName;
		std::ofstream m_LogFile;
	};

}//namespace hamiltonian
}//namespace mcpack

#endif //MCPACK_SAMPLER_HPP
