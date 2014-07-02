#ifndef MCPACK_MPIINTERCHAINCLASSICHMC_HPP
#define MCPACK_MPIINTERCHAINCLASSICHMC_HPP

namespace mcpack{ namespace hamiltonian {

    template<class _ProposalType>
    class MPIInterChainClassicHMC
    {
    public:
        typedef _ProposalType ProposalType;
        typedef typename ProposalType::RealType RealType;
        typedef typename ProposalType::IndexType IndexType;
        typedef typename ProposalType::RealVectorType RealVectorType;
        typedef typename ProposalType::RandVarGenType RandVarGenType;
        typedef typename Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic> RealMatrixType;
        typedef typename RandVarGenType::SeedType SeedType;

        typedef boost::mpi::environment MpiEnvType;
        typedef boost::mpi::communicator MpiCommType;
        typedef boost::mpi::request MpiRequestType;

        MPIInterChainClassicHMC()
        :m_q0(0,0),m_AccRate(0),m_RVGen(0)
        {}

        MPIInterChainClassicHMC(ProposalType const & Prop,RealVectorType const& q0,const SeedType seed)
        :m_Prop(Prop),m_q0(q0),m_AccRate(0),m_RVGen(seed)
        {
            GenerateRandomSeed();
        }

        void Generate(RealMatrixType & Samples)
        {
            IndexType NSamples=Samples.rows();
            IndexType NDim=Samples.cols();

            //MCPACK_ASSERT(m_Discr.NDim()==NDim,
                //"DiscretisationType and Samples should have the same dimensionality.");

            IndexType iter=0;
            IndexType samp=0;

            while(samp < NSamples)
            {
                //if(m_World.rank() == 0 ) std::cout<<"sample "<<samp<<std::endl;

                RealVectorType q1(m_q0);
                RealVectorType DeltaH=RealVectorType::Zero(m_World.size());

                //propose. ie. each process proposes a new point
                //if(m_World.rank() == 0 ) std::cout<<"making proposals "<<std::endl;
                for(IndexType i=0;i<m_World.size();++i)
                {
                    if(i == (IndexType) m_World.rank())
                    {
                        RealType dH=0;
                        m_Prop.Propose(q1,dH,m_RVGen);
                        DeltaH(i)=dH;
                    }
                }

                
                m_World.barrier();
                //if(m_World.rank() == 0 ) std::cout<<"done making proposals now distributing"<<std::endl;


                //distribute the computed deltaH
                DistributeDeltaH(DeltaH);

                m_World.barrier();


                //if(m_World.rank() == 0 ) std::cout<<"done distributing now finding min(deltaH)"<<std::endl;

                //now find the minimum of deltaH vector
                IndexType dhInd=FindMinDeltaH(DeltaH);
                RealType dHMin=DeltaH(dhInd);

                //std::cout<<"\n process "<<m_World.rank()<<"\n"<<DeltaH<<"\n"<<std::endl;
                m_World.barrier();

                //std::cout<<"\n process "<<m_World.rank()<<"\n"<<dhInd<<"\n"<<std::endl;

                
                //if(m_World.rank() == 0 ) std::cout<<"done finding min(deltaH) now accept/reject"<<std::endl;

                RealType u=m_RVGen.Uniform();

                DistributeProposal(q1,0);

                m_World.barrier();
                std::cout<<"process "<<m_World.rank()<<" "<<dHMin<<std::endl;
                m_World.barrier();
                if(u < exp(-dHMin))
                {
                    //std::cout<<"process "<<m_World.rank()<<" accepted"<<std::endl;
                    //distribute the accepted proposal 
                    //to all other processes
                    //if(m_World.rank() == 0 ) std::cout<<"done accept/reject now distributing proposal dhInd= "<<dhInd<<std::endl;                
                    
                    //DistributeProposal(q1,0);

                    //if(m_World.rank() == 0 ) std::cout<<"done distributing proposal"<<std::endl;  

                    
                    //m_q0=q1;
                    //Samples.row(samp)=m_q0;
                    //++samp;
                    samp=samp+1;
                    m_World.barrier();
                    double x=10;
                    
                }
                
                
                
                iter++;

                if(iter > 1000 ) break;

                m_World.barrier();
                //std::cout<<std::endl;
            }

            //m_AccRate=(RealType)samp/(RealType)iter;
        }

        void DistributeDeltaH(RealVectorType & DeltaH) const
        {

            /*
            IndexType rind=0;
            std::vector<MpiRequestType> reqs1(m_World.size()-1);

            for(IndexType source=0;source<  (IndexType)  m_World.size();++source)
            {                    
                if(source ==  (IndexType) m_World.rank() )
                {
                    // send my dh to every other process
                    std::vector<MpiRequestType> reqs2(m_World.size()-1);

                    IndexType sind=0;
                    for(IndexType dest=0;dest<  (IndexType)  m_World.size();++dest)
                    {
                        if(dest != (IndexType) m_World.rank() )
                        {
                            reqs2[sind] =  m_World.isend(dest, source, DeltaH(source));
                            ++sind;
                        }                            
                    }

                    boost::mpi::wait_all(reqs2.begin(), reqs2.end());
                }
                else
                {
                    //receive the values of all ther dhs from other processes
                    reqs1[rind] =  m_World.irecv(source, source, DeltaH(source));
                    ++rind;
                }
            }
            boost::mpi::wait_all(reqs1.begin(), reqs1.end());

            */

            //send my dh to every other process
            IndexType source=m_World.rank();

            std::vector<MpiRequestType> reqs2(m_World.size()-1);
            IndexType sind=0;

            //if(m_World.rank() == 0 ) std::cout<<"@DistributeDeltaH sending to destimations"<<std::endl;  
            for(IndexType dest=0;dest<  (IndexType)  m_World.size();++dest)
            {                
                if(dest != (IndexType) m_World.rank() ) //I shouldn't be the destimation myself
                {
                    //reqs2[sind] =  m_World.isend(dest, source, DeltaH(source));//send the current dh to every other process
                    m_World.send(dest, source, DeltaH(source));//send the current dh to every other process
                    ++sind;
                }
            }
            //boost::mpi::wait_all(reqs2.begin(), reqs2.end());

            //m_World.barrier();

            //if(m_World.rank() == 0 ) std::cout<<"@DistributeDeltaH receiving from sources"<<std::endl;
            //receive dh from every other process
            std::vector<MpiRequestType> reqs1(m_World.size()-1);
            IndexType rind=0;
            for(source=0;source<  (IndexType)  m_World.size();++source)
            {
                if(source !=  (IndexType) m_World.rank() )//I shouldn't reveive stuff from myself
                {
                    //reqs1[rind] =  m_World.irecv(source, source, DeltaH(source));//get the values for the indices other than my index
                    m_World.recv(source, source, DeltaH(source));//get the values for the indices other than my index
                    ++rind;
                }
            }
            //boost::mpi::wait_all(reqs1.begin(), reqs1.end());
        }

        IndexType FindMinDeltaH(RealVectorType const& DeltaH) const
        {
            RealType minDh=DeltaH(0);
            IndexType minInd=0;
            for(IndexType i=1;i<DeltaH.rows();++i)
            {
                if(DeltaH(i) < minDh)
                {
                    minDh = DeltaH(i);
                    minInd = i;
                }
            }

            return minInd;
        }

        void DistributeProposal(RealVectorType & q1,const IndexType source) const
        {
            /*
            if((IndexType) m_World.rank() == source)
            {
                std::vector<MpiRequestType> reqs1(m_World.size()-1);
                IndexType dind=0;
                for(IndexType  dest=0;dest< (IndexType) m_World.size();++dest)
                {
                    if(dest != source )
                    {
                        reqs1[dind] =  m_World.isend(dest, source, &q1(0,0),q1.rows());
                        ++dind;
                    }
                }

                boost::mpi::wait_all(reqs1.begin(), reqs1.end());
            }
            else
            {
                std::vector<MpiRequestType> reqs2(m_World.size()-1);
                IndexType dind = m_World.rank() > source ? m_World.rank()-1 : m_World.rank();
                reqs2[ dind ] =  m_World.irecv(source, source, &q1(0,0),q1.rows());
                boost::mpi::wait_all(reqs2.begin(), reqs2.end());
            }
            */

            if((IndexType) m_World.rank() == source)
            {
                //std::vector<MpiRequestType> reqs1(m_World.size()-1);
                IndexType dind=0;
                for(IndexType  dest=0;dest< (IndexType) m_World.size();++dest)
                {
                    if(dest != source )
                    {
                        //double x=0;
                        //reqs1[dind] =  m_World.isend(dest, source, &q1(0,0),q1.rows());
                        //reqs1[dind] =  m_World.isend(dest, source, x);
                        m_World.send(dest, source, &q1(0,0),q1.rows());
                        ++dind;
                    }
                }
                //boost::mpi::wait_all(reqs1.begin(), reqs1.end());
            }
            else
            {
                //MpiRequestType req = m_World.irecv(source, source, &q1(0,0),q1.rows());
                //double x=0;
                //MpiRequestType req = m_World.irecv(source, source, x);
                m_World.recv(source, source, &q1(0,0),q1.rows());
                //boost::mpi::wait_all(&req,&req);
            }

            //m_World.barrier();
        }

        RealType GetAcceptanceRate(void) const
        {
            return m_AccRate;
        }

        //the folling methods are required for 
        //setting seeds from outside
        void SetSeed(unsigned long seed)
        {
            m_RVGen.Seed(seed);
        }

        void GetRandState(std::stringstream & RNGstate) const
        {
            m_RVGen.GetState(RNGstate);
        }

        void SetRandState(std::stringstream & RNGstate)
        {
            m_RVGen.SetState(RNGstate);
        }

        void SetStartPoint(RealVectorType const & q0)
        {
            m_q0=q0;
        }

        void GenerateRandomSeed()
        {
            std::vector<unsigned long> seedVect(m_World.size());
            for(size_t i=0;i<(size_t)m_World.size();++i)
            {
                seedVect[i]=(unsigned long)rand();
            }
            m_RVGen.Seed(seedVect[m_World.rank()]);
        }

    private:
        //PropVectorType m_PropVector;
        //DeltaHVectorType m_DeltaVector;
        ProposalType m_Prop;
        RealVectorType m_q0;
        RealType m_AccRate;
        RandVarGenType m_RVGen;
        MpiEnvType m_Env;
        MpiCommType m_World;
    };

}//namespace hamiltonian
}//namespace mcpack

#endif //MCPACK_MPIINTERCHAINCLASSICHMC_HPP