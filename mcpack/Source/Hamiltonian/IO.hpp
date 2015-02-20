/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_IO_HPP
#define MCPACK_IO_HPP

namespace mcpack { namespace hamiltonian{

    /**
     * \ingroup Hamiltonian
     *
     * \class IOWriteAllParams
     *
     * \brief A class for managing the IO.
     *
     * \tparam _realMatrixType real matrix type
     *
     * This class manages the IO of the sampler. In particular this class
     * writes the states of all the parameters to the chain output file.
     */
    template<class _realMatrixType>
    class IOWriteAllParams
    {
    public:

        /**
         * \typedef _realMatrixType realMatrixType
         * \brief real matrix type
         */
        typedef _realMatrixType realMatrixType;

        /**
         * \typedef typename realMatrixType::Scalar realScalarType
         * \breif floating point type
         */
        typedef typename realMatrixType::Scalar realScalarType;

        /**
         * \typedef typename realMatrixType::Index indexType
         * \breif integral type
         */
        typedef typename realMatrixType::Index indexType;

        static_assert(std::is_floating_point<realScalarType>::value,
            "PARAMETER SHOULD BE A FLOATING POINT TYPE");

        /**
         * \brief The default constructor
         */
        IOWriteAllParams(){}

        /**
         * \brief A constructor that sets up the IO for the sampler
         * \param  fileName The name of the output file
         */
        explicit IOWriteAllParams(std::string const& fileName)
        :m_fileName(fileName),m_delimiter(std::string(",")),m_precision(10)
        {

        }

        /**
         * \brief A constructor that sets up the IO for the sampler
         * \param fileName The name of the output file
         * \param delimiter The delimiter between two numbers e.g. a coma
         * \param precision The precision with which the output is written
         */
        IOWriteAllParams(std::string const& fileName,std::string const & delimiter,
            unsigned int const precision)
        :m_fileName(fileName),m_delimiter(delimiter),m_precision(precision)
        {

        }

        /**
         * \brief The copy constructor
         * \param other The class from which the members are to be copied
         */
        IOWriteAllParams(IOWriteAllParams  const & other)
        {
            m_fileName = other.m_fileName;
            m_delimiter = other.m_delimiter;
            m_precision = other.m_precision;
        }

        /**
         * \breif The default destructor
         */
        ~IOWriteAllParams()
        {
            m_file.close();
        }

        /**
         * \brief write the samples to a file
         * \param samples the MCMC samples from the sampler
         */
        void write(realMatrixType const& samples)
        {
            if(!m_file.is_open())
            {
                m_file.open(m_fileName.c_str(),std::ios::app);
            }

            if(m_file.is_open())
            {
                m_file<<std::scientific;
                m_file<<std::setprecision(m_precision);
                for(indexType i=0;i<samples.rows();++i)
                {
                    for(indexType j=0;j<samples.cols()-1;++j)
                    {
                        m_file<<samples(i,j)<<m_delimiter;
                    }
                    m_file<<samples(i,(samples.cols()-1) )<<std::endl;
                }
            }
            else
            {
                std::string message = std::string("Error in opening the file ")+m_fileName;
                throw std::runtime_error(message);
            }
        }

        /**
         * \brief Return the output file name
         * \return  the output file name
         */
        inline std::string getFileName(void) const
        {
            return m_fileName;
        }

        /**
         * \brief Set the output file name
         * \param fileName the output file name
         */
        inline void setFileName(std::string const& fileName)
        {
            m_fileName = fileName;
        }

    private:
        std::string m_fileName; /**< the file output file name */
        std::ofstream m_file; /**< output file handler*/
        std::string m_delimiter; /**< delimiter */
        unsigned int m_precision; /**< precision */

    };
}
}

#endif //MCPACK_IO_HPP
