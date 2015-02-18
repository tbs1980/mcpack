/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_RUNTIMECONTROL_HPP
#define MCPACK_RUNTIMECONTROL_HPP


namespace mcpack { namespace hamiltonian {


    /**
     * \ingroup Hamiltonian
     *
     * \class runCtrlFiniteSamples
     *
     * \brief A class for controlling runtime behaviour of the sampler
     *
     * \tparam _realMatrixType real matrix type
     *
     * This class controls the runtime behaviour of the sampler.
     * TODO more info here.
     */
    template<class _realMatrixType>
    class runCtrlFiniteSamples
    {
    public:

        /**
         * \typedef _realMatrixType realMatrixType
         * \brief real matrix type
         */
        typedef _realMatrixType realMatrixType;

        /**
         * \typedef typename realMatrixType::Index indexType
         * \breif integral type
         */
        typedef typename realMatrixType::Index indexType;

        /**
         * \typedef typename realMatrixType::Scalar realScalarType
         * \brief floating point type
         */
        typedef typename realMatrixType::Scalar realScalarType;

        static_assert(std::is_floating_point<realScalarType>::value,
            "PARAMETER SHOULD BE A FLOATING POINT TYPE");

        /**
         * \brief The default constructor
         */
        runCtrlFiniteSamples()
        :m_numParams(0),m_samples(0),m_numSamples(0),
        m_packetSize(0),m_burn(0),m_numBurn(0)
        {

        }

        /**
         * \brief A constructor that sets up the runtime control
         * \param numParams The number of parameters
         * \param numSamples The number of samples required
         * \param packetSize The size of the packet of samples for each iteration
         * \param numBurn The number of samples to be burned / ignored
         * \param root The root of the path to output files
         * \param silent Output to console?
         */
        runCtrlFiniteSamples(indexType const numParams,indexType const numSamples,
            indexType const packetSize, indexType const numBurn,
            std::string const& root, bool silent)
        :m_numParams(numParams),m_samples(0),m_numSamples(numSamples),
        m_packetSize(packetSize),m_burn(0),m_numBurn(numBurn),m_root(root),
        m_silent(silent),m_logFileName(root+std::string(".log")),m_resume(false),
        m_continueSampling(true),m_logFileHasHeader(false)
        {
            BOOST_ASSERT_MSG(m_numSamples>0,"Maximum number of samples should be a positive integer");
            BOOST_ASSERT_MSG(m_numBurn>=0,"Number of samples to be burned should be a >= 0");
        }

        /**
         * \brief A function to save the control parameters
         * \param samples   The matrix of samples
         * \param randState The random number state
         * \param accRate   acceptance rate
         */
        void save(realMatrixType const & samples,std::stringstream const& randState,
            realScalarType const accRate)
        {
            m_randState = randState.str();

            indexType n = samples.rows();
            if(m_burn >= m_numBurn)
            {
                m_samples += n;
            }
            else
            {
                m_burn += n;
            }

            m_continueSampling = m_samples >= m_numSamples ? false : true;

            if(!m_logFileHasHeader)
            {
                writeInfo2LogFile();
            }

            m_pt.put("control.burn",(indexType) m_burn);
            m_pt.put("control.samples",(indexType) m_samples);
            m_pt.put("chain.accRate",(realScalarType) accRate);

            std::stringstream chainState;
            for(indexType i=0;i<samples.cols()-1;++i)
            {
                chainState<<samples(samples.rows()-1,i)<<" ";
            }
            chainState<<samples(samples.rows()-1,samples.cols()-1);
            m_chainState = chainState.str();
            m_pt.put("chain.State",(std::string) chainState.str());
            m_pt.put("random.State",(std::string) randState.str());

            // TODO make this fine and xml file instead?
            boost::property_tree::ini_parser::write_ini(m_logFileName,m_pt);
        }

        /**
         * \brief Return the number of parameters
         * \return  the number of parameters
         */
        inline indexType numParams(void) const
        {
            return m_numParams;
        }

        /**
         * \brief Return the packet size
         * \return  the packet size
         */
        inline indexType packetSize(void) const
        {
            return m_packetSize;
        }

        /**
         * \brief Return the root path to the output files
         * \return  the root path to the output files
         */
        inline std::string const& root(void) const
        {
            return m_root;
        }

        /**
         * \brief Return the logical paramter for console output
         * @return  the logical paramter for console output
         */
        inline bool silent(void) const
        {
            return m_silent;
        }

        /**
         * \brief Return the logical parameter for resume
         * \return  the logical parameter for resume
         */
        inline bool resume(void) const
        {
            return m_resume;
        }

        /**
         * \breif Return the random number state
         * @return  the random number state
         */
        inline std::string const & randState(void) const
        {
            return m_randState;
        }

        /**
         * \breif Return the chain state
         * @return  the chain state
         */
        inline std::string const & chainState(void) const
        {
            return m_chainState;
        }

        /**
         * \breif Return the log file name
         * @return  the log file name
         */
        inline std::string const & getLogFileName(void) const
        {
            return m_logFileName;
        }

        /**
         * \brief Set the log file name
         * \param LogFileName the log file name
         */
        inline void setLogFileName(std::string const & LogFileName)
        {
            m_logFileName=LogFileName;
        }

        /**
         * \breif Return the logical parameter for continuing sampling
         * @return the logical parameter for continuing sampling
         */
        inline bool continueSampling(void) const
        {
            return m_continueSampling;
        }

        /**
         * \brief Load the information from the log file
         */
        void loadInfoFromLogFile(void)
        {
            try
            {
                boost::property_tree::ini_parser::read_ini(m_logFileName,m_pt);

                //now assign data
                m_numParams=m_pt.get<indexType>("control.numParams");
                m_samples=m_pt.get<indexType>("control.samples");
                m_numSamples=m_pt.get<indexType>("control.numSamples");
                m_packetSize=m_pt.get<indexType>("control.packetSize");
                m_burn=m_pt.get<indexType>("control.burn");
                m_numBurn=m_pt.get<indexType>("control.numBurn");
                m_root=m_pt.get<std::string>("control.root");
                m_silent=m_pt.get<bool>("control.Silent");
                m_chainState=m_pt.get<std::string>("chain.State");
                m_randState=m_pt.get<std::string>("random.State");
                m_resume=true;
                m_continueSampling = m_samples >= m_numSamples ? false : true;

            }
            catch(std::exception& e)
            {
                m_resume=false;
            }

        }

        /**
         * \brief Write the information to a log file
         */
        void writeInfo2LogFile(void)
        {
            m_pt.put("control.numParams",(indexType) m_numParams);
            m_pt.put("control.numSamples",(indexType) m_numSamples);
            m_pt.put("control.numBurn",(indexType) m_numBurn);
            m_pt.put("control.packetSize",(indexType) m_packetSize);
            m_pt.put("control.root",(std::string)  m_root);
            m_pt.put("control.Silent",(indexType) m_silent);

            boost::property_tree::ini_parser::write_ini(m_logFileName,m_pt);

            m_logFileHasHeader = true;
        }

    private:
        indexType m_numParams; /**< number of parameters */
        indexType m_samples; /**< number of samples taken so far*/
        indexType m_numSamples; /**< number samples required */
        indexType m_packetSize; /**< packet size */
        indexType m_burn; /**< number of samples burned so far*/
        indexType m_numBurn; /**< number samples to be burned */
        std::string m_root; /**< root path to log files */
        bool m_silent; /**< output to console */
        std::string m_logFileName; /**< log file name */
        boost::property_tree::ptree m_pt; /**< a tree for control parameters */
        std::string m_randState; /**< random number state */
        std::string m_chainState; /**< chain state */
        bool m_resume; /**< resume from previous state */
        bool m_continueSampling; /**< continue sampling */
        bool m_logFileHasHeader; /**< log file header present */
    };

}
}

#endif //MCPACK_RUNTIMECONTROL_HPP
