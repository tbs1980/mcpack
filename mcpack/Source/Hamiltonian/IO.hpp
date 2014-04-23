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

#ifndef MCPACK_IO_HPP
#define MCPACK_IO_HPP

namespace mcpack { namespace hamiltonian{

	template<class _MatrixType>
	class IO_WriteAll
	{
	public:
		typedef _MatrixType MatrixType;
		typedef typename MatrixType::Index IndexType;

		explicit IO_WriteAll(std::string const FileName)
		:m_FileName(FileName),m_Separation(std::string(",")),m_precision(10)
		{
			m_File.open(m_FileName.c_str(),std::ios::trunc);

			if(!m_File.is_open())
			{
				std::string message=std::string("Error in opening the file ")+m_FileName;
				throw mcpack::utils::TextDataException(message);			
			}

		}

		~IO_WriteAll()
		{
			m_File.close();
		}
	
		void Write(MatrixType const& Samples)
		{
			if(m_File.is_open())
			{
				m_File<<std::scientific;
				m_File<<std::setprecision(m_precision);
				for(IndexType i=0;i<Samples.rows();++i)
				{
					for(IndexType j=0;j<Samples.cols()-1;++j)
					{
						m_File<<Samples(i,j)<<m_Separation;
					}
					m_File<<Samples(i,(Samples.cols()-1) )<<std::endl;
				}				
			}
			else
			{
				std::string message=std::string("Error in opening the file ")+m_FileName;
				throw mcpack::utils::TextDataException(message);
			}
		}
	private:
		std::string m_FileName;
		std::ofstream m_File;
		std::string m_Separation;
		unsigned m_precision;
		
	};
}
}

#endif //MCPACK_IO_HPP
