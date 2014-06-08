/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_IO_HPP
#define MCPACK_IO_HPP

namespace mcpack { namespace hamiltonian{

    template<class _MatrixType>
    class IO_WriteAll
    {
    public:
        typedef _MatrixType MatrixType;
        typedef typename MatrixType::Index IndexType;

        IO_WriteAll(){}

        explicit IO_WriteAll(std::string const FileName)
        :m_FileName(FileName),m_Separation(std::string(",")),m_precision(10)
        {

        }

        IO_WriteAll(IO_WriteAll  const & other)
        {
            m_FileName=other.m_FileName;
            m_Separation=other.m_Separation;
            m_precision=other.m_precision;
        }

        ~IO_WriteAll()
        {
            m_File.close();
        }
    
        void Write(MatrixType const& Samples)
        {
            if(!m_File.is_open())
            {
                m_File.open(m_FileName.c_str(),std::ios::app);
            }

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

        std::string GetFileName(void) const
        {
            return m_FileName;
        }

        void SetFileName(std::string const& FileName)
        {
            m_FileName=FileName;
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
