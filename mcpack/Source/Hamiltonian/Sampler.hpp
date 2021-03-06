/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
        typedef typename EngineType::RealType RealType;
        typedef typename EngineType::IndexType IndexType;
        
        Sampler(EngineType const & Eng,IOType const& IO,RunCtrlType const& RunCtrl)
        :m_Eng(Eng),m_IO(IO),m_RunCtrl(RunCtrl)
        {
            m_RunCtrl.LoadInfoFromLogFile();

            if(!m_RunCtrl.Resume())
            {
                WriteOutput2Console(std::string("No resume files present "),m_RunCtrl.Silent());

                m_RunCtrl.WriteInfo2LogFile();
            }
        }

        void Run()
        {
            RealMatrixType Samples(m_RunCtrl.PacketSize(),m_RunCtrl.NumParas());

            if(m_RunCtrl.Resume())
            {
                WriteOutput2Console(std::string("Resuming from previous run"),m_RunCtrl.Silent());

                if(m_RunCtrl.Continue())
                {
                    std::stringstream RandState;
                    RandState<<m_RunCtrl.RandState();
                    m_Eng.SetRandState(RandState);

                    RealVectorType ChainState=
                        mcpack::utils::String2Vector<RealVectorType>(m_RunCtrl.ChainState(),std::string(" "));
                    m_Eng.SetStartPoint(ChainState);
                }   
                else
                {
                    WriteOutput2Console(std::string("Sampling already finished in the previous run"),m_RunCtrl.Silent());
                }
            }
            else
            {
                WriteOutput2Console(std::string("Starting sampling from scratch"),m_RunCtrl.Silent());
            }
            
            while(m_RunCtrl.Continue())
            {
                std::stringstream RandState;

                m_Eng.Generate(Samples);

                m_Eng.GetRandState(RandState);
                RealType AccRate=m_Eng.GetAcceptanceRate();

                m_RunCtrl.Save(Samples,RandState,AccRate);

                m_IO.Write(Samples);
            }
        }

        static void WriteOutput2Console(std::string message,bool silent)
        {
            if(!silent)
            {
                std::cout<<"--> "<<message<<"\n"<<std::endl;
            }
        }

    private:

        EngineType m_Eng;
        IOType m_IO;
        RunCtrlType m_RunCtrl;
    };

}//namespace hamiltonian
}//namespace mcpack

#endif //MCPACK_SAMPLER_HPP
