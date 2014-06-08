/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#ifndef MCPACK_MPISAMPLER_HPP
#define MCPACK_MPISAMPLER_HPP

namespace mcpack{ namespace hamiltonian {

    template<class _Engine,class _IOType,class _RunCtrlType>
    class Mpi_Sampler
    {
    public:
        typedef _Engine EngineType;
        typedef _IOType IOType;
        typedef _RunCtrlType RunCtrlType;

        typedef typename EngineType::RealMatrixType RealMatrixType;
        typedef typename EngineType::RealVectorType RealVectorType;
        typedef typename EngineType::RealType RealType;
        typedef typename EngineType::IndexType IndexType;

        typedef boost::mpi::environment MpiEnvType;
        typedef boost::mpi::communicator MpiCommType;
        typedef boost::mpi::request MpiRequestType;
        
        Mpi_Sampler(EngineType const & Eng,IOType const& IO,RunCtrlType const& RunCtrl)
        :m_Eng(Eng),m_IO(IO),m_RunCtrl(RunCtrl)
        {
            std::stringstream ss;
            ss<<"."<<m_World.rank();

            std::string IOFileName=m_IO.GetFileName();
            IOFileName+=ss.str();
            m_IO.SetFileName(IOFileName);

            std::string LogFileName=m_RunCtrl.GetLogFileName();
            LogFileName+=ss.str();
            m_RunCtrl.SetLogFileName(LogFileName);

            m_RunCtrl.LoadInfoFromLogFile();

            if(!m_RunCtrl.Resume())
            {
                m_RunCtrl.WriteInfo2LogFile();
            }

            PrintInfo();
        }

        void Run()
        {
            RealMatrixType Samples(m_RunCtrl.PacketSize(),m_RunCtrl.NumParas());
            
            if(m_RunCtrl.Resume() )
            {
                if(m_RunCtrl.Continue())
                {
                    std::stringstream RandState;
                    RandState<<m_RunCtrl.RandState();
                    m_Eng.SetRandState(RandState);

                    RealVectorType ChainState=
                        mcpack::utils::String2Vector<RealVectorType>(m_RunCtrl.ChainState(),std::string(" "));

                    m_Eng.SetStartPoint(ChainState);                    
                }
            }
            else
            {
                GenerateRandomSeed();               
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


    private:


        void GenerateRandomSeed()
        {
            std::vector<unsigned long> seedVect(m_World.size());
            for(size_t i=0;i<(size_t)m_World.size();++i)
            {
                seedVect[i]=(unsigned long)rand();
            }
            m_Eng.SetSeed(seedVect[m_World.rank()]);
        }

        void PrintInfo(void) const
        {

            //print a summary of what we have found
            if (m_World.rank() == 0)
            {
                std::string Header("\n----------------------------------------------------------------\n"); 
                Header+=std::string("        Chain        Logfile present?        Sampling finished?\n");
                Header+=std::string("----------------------------------------------------------------");

                WriteOutput2Console(Header);

                std::vector<MpiRequestType> reqs(m_World.size()-1);

                //Is chain zero resuming from previous run?
                std::string IsResuming=Bool2String(m_RunCtrl.Resume());
                std::string IsSampling=Bool2String(!m_RunCtrl.Continue());

                std::stringstream ss_resume;
                ss_resume<<"         0                 "<<IsResuming<<"                     "<<IsSampling;
                WriteOutput2Console(ss_resume.str());

                //we need messages from other np-1 processes.
                std::vector<std::string> msgVect(m_World.size()-1);

                //receive the messages from other processes
                for(IndexType i=1;i<m_World.size();++i)
                {
                    reqs[i-1] =  m_World.irecv(i, i, msgVect[i-1]);
                }

                boost::mpi::wait_all(reqs.begin(), reqs.end());

                for(IndexType i=1;i<m_World.size();++i)
                {
                    std::stringstream SSResumeOthers;
                    SSResumeOthers<<"         "<<i<<"                 "<<msgVect[i-1];
                    WriteOutput2Console(SSResumeOthers.str());
                }
            }
            else
            {
                std::vector<MpiRequestType> reqs(m_World.size()-1);

                std::string IsResuming=Bool2String(m_RunCtrl.Resume());
                std::string IsSampling=Bool2String(!m_RunCtrl.Continue());
                std::string messg=IsResuming+std::string("                     ")+IsSampling;

                //send the messages to process 0
                reqs[m_World.rank()-1] =  m_World.isend(0, m_World.rank(), messg);

                boost::mpi::wait_all(reqs.begin(), reqs.end());
            }
        }

        std::string Bool2String(bool val) const
        {
            if(val==true)
            {
                return std::string("Yes");
            }
            else
            {
                return std::string("No");
            }
        }
        
        void WriteOutput2Console(std::string message) const
        {
            if(!m_RunCtrl.Silent())
            {
                std::cout<<message<<std::endl;
            }
        }

        EngineType m_Eng;
        IOType m_IO;
        RunCtrlType m_RunCtrl;
        MpiEnvType m_Env;
        MpiCommType m_World;

    };

}//namespace hamiltonian
}//namespace mcpack

#endif //MCPACK_MPISAMPLER_HPP