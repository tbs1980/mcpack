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
		typedef typename MatrixType::Scalar RealType;

		RunCtrl_FiniteSamples()
		:m_NumParas(0),m_Samples(0),m_NumSamples(0),
		m_PacketSize(0),m_Burn(0),m_NumBurn(0)
		{}

		RunCtrl_FiniteSamples(IndexType const NumParas,IndexType const NumSamples, 
			IndexType const PacketSize,	IndexType const NumBurn,std::string const& root,
			bool silent)
		:m_NumParas(NumParas),m_Samples(0),m_NumSamples(NumSamples),
		m_PacketSize(PacketSize),m_Burn(0),m_NumBurn(NumBurn),m_Root(root),
		m_Silent(silent),m_LogFileName(root+std::string(".log")),m_Resume(false),
		m_Continue(true),m_LogFileHasHeader(false)
		{
			MCPACK_ASSERT(m_NumSamples>0,"Maximum number of samples should be a positive integer");
			MCPACK_ASSERT(m_NumBurn>=0,"Number of samples to be burned should be a >= 0");
		}

		void Save(MatrixType const & Samples,std::stringstream const& RandState,
			RealType const AccRate)
		{
			m_RandState=RandState.str();

			IndexType n=Samples.rows();
			if(m_Burn >= m_NumBurn)
			{
				m_Samples+=n;
			}
			else
			{
				m_Burn+=n;
			}

			m_Continue = m_Samples >= m_NumSamples ? false : true;

			if(!m_LogFileHasHeader)
			{
				WriteInfo2LogFile();
			}

			m_Pt.put("Control.Burn",(IndexType) m_Burn);
			m_Pt.put("Control.Samples",(IndexType) m_Samples);
			m_Pt.put("Chain.AccRate",(RealType) AccRate);

			std::stringstream ChainState;
			for(IndexType i=0;i<Samples.cols()-1;++i)
			{
				ChainState<<Samples(Samples.rows()-1,i)<<" ";
			}
			ChainState<<Samples(Samples.rows()-1,Samples.cols()-1);
			m_ChainState=ChainState.str();
			m_Pt.put("Chain.State",(std::string) ChainState.str());
			m_Pt.put("Random.State",(std::string) RandState.str());

			boost::property_tree::ini_parser::write_ini(m_LogFileName,m_Pt);
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
			return m_Root;
		}

		bool Silent(void) const
		{
			return m_Silent;
		}

		bool Resume(void) const
		{
			return m_Resume;
		}

		std::string const & RandState(void) const
		{
			return m_RandState;
		}

		std::string const & ChainState(void) const
		{
			return m_ChainState;
		}

		std::string const & GetLogFileName(void) const
		{
			return m_LogFileName;
		}

		void SetLogFileName(std::string const & LogFileName)
		{
			m_LogFileName=LogFileName;
		}

		bool Continue() const
		{
			return m_Continue;
		}

		void LoadInfoFromLogFile()
		{
			try
			{
				boost::property_tree::ini_parser::read_ini(m_LogFileName,m_Pt);

				//now assign data
				m_NumParas=m_Pt.get<IndexType>("Control.NumParas");
				m_Samples=m_Pt.get<IndexType>("Control.Samples");
				m_NumSamples=m_Pt.get<IndexType>("Control.NumSamples");
				m_PacketSize=m_Pt.get<IndexType>("Control.PacketSize");
				m_Burn=m_Pt.get<IndexType>("Control.Burn");
				m_NumBurn=m_Pt.get<IndexType>("Control.NumBurn");
				m_Root=m_Pt.get<std::string>("Control.Root");
				m_Silent=m_Pt.get<bool>("Control.Silent");
				m_ChainState=m_Pt.get<std::string>("Chain.State");
				m_RandState=m_Pt.get<std::string>("Random.State");
				m_Resume=true;
				m_Continue = m_Samples >= m_NumSamples ? false : true;

			}
			catch(std::exception& e)
			{
				m_Resume=false;
			}
			
		}

		void WriteInfo2LogFile(void)
		{
			m_Pt.put("Control.NumParas",(IndexType) m_NumParas);
			m_Pt.put("Control.NumSamples",(IndexType) m_NumSamples);
			m_Pt.put("Control.NumBurn",(IndexType) m_NumBurn);
			m_Pt.put("Control.PacketSize",(IndexType) m_PacketSize);
			m_Pt.put("Control.Root",(std::string)  m_Root);
			m_Pt.put("Control.Silent",(IndexType) m_Silent);

			boost::property_tree::ini_parser::write_ini(m_LogFileName,m_Pt);

			m_LogFileHasHeader=true;		
		}

	private:


		IndexType m_NumParas;
		IndexType m_Samples;
		IndexType m_NumSamples;
		IndexType m_PacketSize;
		IndexType m_Burn;
		IndexType m_NumBurn;
		std::string m_Root;
		bool m_Silent;
		std::string m_LogFileName;
		boost::property_tree::ptree m_Pt;
		std::string m_RandState;
		std::string m_ChainState;
		bool m_Resume;
		bool m_Continue;
		bool m_LogFileHasHeader;
	};

}
}

#endif //MCPACK_RUNTIMECONTROL_HPP