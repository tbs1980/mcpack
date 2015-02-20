/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_sampler_HPP
#define MCPACK_sampler_HPP

namespace mcpack{ namespace hamiltonian{

    /**
     * \ingroup Hamiltonian
     *
     * \class sampler
     *
     * \brief A class for performing MCMC sampling
     *
     * \tparam _engineType MCMC Engine type
     * \tparam _IOType Input-output type
     * \tparam _runCtrlType Runtime-control type
     *
     * A class that performs the MCMC sampling using the specified engine,
     * IO and Runtime-control.
     */
    template<class _engineType,class _IOType,class _runCtrlType>
    class sampler
    {
    public:

        /**
         * \typedef _engineType engineType
         * \brief MCMC Engine type
         */
        typedef _engineType engineType;

        /**
         * \typedef _IOType IOType
         * \brief Input-output type
         */
        typedef _IOType IOType;

        /**
         * \typedef _runCtrlType runCtrlType
         * \brief Runtime-control type
         */
        typedef _runCtrlType runCtrlType;

        /**
         * \typedef typename engineType::realMatrixType realMatrixType
         * \brief real matrix type
         */
        typedef typename engineType::realMatrixType realMatrixType;

        /**
         * \typedef typename engineType::realVectorType realVectorType
         * \breif real vector type
         */
        typedef typename engineType::realVectorType realVectorType;

        /**
         * \typedef typename engineType::realScalarType realScalarType
         * \brief real floating point type
         */
        typedef typename engineType::realScalarType realScalarType;

        /**
         * \typedef typename engineType::indexType indexType
         * \brief integral type
         */
        typedef typename engineType::indexType indexType;

        static_assert(std::is_floating_point<realScalarType>::value,
            "PARAMETER SHOULD BE A FLOATING POINT TYPE");

        /**
         * \brief A constructor that sets up the sampler
         *
         * \param eng MCMC Engine
         * \param IO input output handle
         * \param runCtrl Runtime-control
         */
        sampler(engineType const & eng,IOType const& IO,runCtrlType const& runCtrl)
        :m_eng(eng),m_IO(IO),m_runCtrl(runCtrl)
        {
            m_runCtrl.loadInfoFromLogFile();

            //TODO make the console output using boost.log

            if(!m_runCtrl.resume())
            {
                writeOutput2Console(std::string("No resume files present "),m_runCtrl.silent());

                m_runCtrl.writeInfo2LogFile();
            }
        }

        /**
         * \brief Run the MCMC sampler
         */
        void run(void)
        {
            realMatrixType samples(m_runCtrl.packetSize(),m_runCtrl.numParams());

            if(m_runCtrl.resume())
            {
                writeOutput2Console(std::string("Resuming from previous run"),m_runCtrl.silent());

                if(m_runCtrl.continueSampling())
                {
                    std::stringstream randState;
                    randState<<m_runCtrl.randState();
                    m_eng.setRandState(randState);

                    //TODO change this to a boost tokeniser?
                    realVectorType chainState=
                        mcpack::utils::String2Vector<realVectorType>(m_runCtrl.chainState(),std::string(" "));
                    m_eng.setStartPoint(chainState);
                }
                else
                {
                    writeOutput2Console(std::string("Sampling already finished in the previous run"),m_runCtrl.silent());
                }
            }
            else
            {
                writeOutput2Console(std::string("Starting sampling from scratch"),m_runCtrl.silent());
            }

            while(m_runCtrl.continueSampling())
            {
                std::stringstream randState;

                m_eng.generate(samples);

                m_eng.getRandState(randState);
                realScalarType accRate=m_eng.getAcceptanceRate();

                m_runCtrl.save(samples,randState,accRate);

                m_IO.write(samples);

                //TODO add tuning facility
                //while( m_runCtrl.continueTuning() )
                //{
                //    m_tuner.tune(samples,accRate)
                //}
            }
        }

        /**
         * \brief write output to console
         * \param message message to output
         * \param silent  is working in silent mode?
         */
        static void writeOutput2Console(std::string message,bool silent)
        {
            if(!silent)
            {
                std::cout<<"--> "<<message<<"\n"<<std::endl;
            }
        }

    private:
        engineType m_eng; /**< MCMC Engine*/
        IOType m_IO; /**< Input-Output Handle*/
        runCtrlType m_runCtrl; /**< Runtime-Control */
    };

}//namespace hamiltonian
}//namespace mcpack

#endif //MCPACK_sampler_HPP
