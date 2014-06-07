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

            PrintResumeInfo();
            PrintSamplingInfo();
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

        void PrintResumeInfo(void) const
        {

            //print a summary of what we have found
            if (m_World.rank() == 0)
            {
                std::vector<MpiRequestType> reqs(m_World.size()-1);

                //Is chain zero resuming from previous run?
                std::string IsResuming=Bool2String(m_RunCtrl.Resume());

                std::stringstream ss_resume;
                ss_resume<<"Chain 0 resuming = "<<IsResuming;
                WriteOutput2Console(ss_resume.str(),m_RunCtrl.Silent());

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
                    SSResumeOthers<<"Chain "<<i<<" resuming = "<<msgVect[i-1];
                    WriteOutput2Console(SSResumeOthers.str(),m_RunCtrl.Silent());
                }
            }
            else
            {
                std::vector<MpiRequestType> reqs(m_World.size()-1);

                std::string IsResuming=Bool2String(m_RunCtrl.Resume());

                //send the messages to process 0
                reqs[m_World.rank()-1] =  m_World.isend(0, m_World.rank(), IsResuming);

                boost::mpi::wait_all(reqs.begin(), reqs.end());
            }
        }

        void  PrintSamplingInfo(void) const
        {

            //print a summary of what we have found
            if (m_World.rank() == 0)
            {
                std::vector<MpiRequestType> reqs(m_World.size()-1);

                //Is chain zero resuming from previous run?
                std::string IsResuming=Bool2String(!m_RunCtrl.Continue());

                std::stringstream ss_resume;
                ss_resume<<"Chain 0 sampling finished = "<<IsResuming;
                WriteOutput2Console(ss_resume.str(),m_RunCtrl.Silent());

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
                    SSResumeOthers<<"Chain "<<i<<" sampling finished = "<<msgVect[i-1];
                    WriteOutput2Console(SSResumeOthers.str(),m_RunCtrl.Silent());
                }
            }
            else
            {
                std::vector<MpiRequestType> reqs(m_World.size()-1);

                std::string IsResuming=Bool2String(!m_RunCtrl.Continue());

                //send the messages to process 0
                reqs[m_World.rank()-1] =  m_World.isend(0, m_World.rank(), IsResuming);

                boost::mpi::wait_all(reqs.begin(), reqs.end());
            }
        }

        std::string Bool2String(bool val) const
        {
            if(val==true)
            {
                return std::string("True");
            }
            else
            {
                return std::string("False");
            }
        }
        
        static void WriteOutput2Console(std::string message,bool silent)
        {
            if(!silent)
            {
                std::cout<<"-->"<<message<<"\n"<<std::endl;
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