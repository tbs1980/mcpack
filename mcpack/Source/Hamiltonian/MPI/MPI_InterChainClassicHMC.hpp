#ifndef MCPACK_MPIINTERCHAINCLASSICHMC_HPP
#define MCPACK_MPIINTERCHAINCLASSICHMC_HPP

namespace mcpack{ namespace hamiltonian {

    template<class _ProposalType>
    class ClassicHMC
    {
    public:
        typedef _ProposalType ProposalType;
        typedef typename ProposalType::RealType RealType;
        typedef typename ProposalType::IndexType IndexType;
        typedef typename ProposalType::RealVectorType RealVectorType;
        typedef typename ProposalType::RandVarGenType RandVarGenType;
        typedef typename Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic> RealMatrixType;
        typedef typename RandVarGenType::SeedType SeedType;
        std::vector<ProposalType> PropVectorType;
        std::vector<RealType> DeltaHVectorType;
        typedef boost::mpi::environment MpiEnvType;
        typedef boost::mpi::communicator MpiCommType;

        ClassicHMC()
        :m_q0(0,0),m_AccRate(0),m_RVGen(0)
        {}

        ClassicHMC(ProposalType const & Prop,RealVectorType const& q0,const SeedType seed)
        :m_q0(q0),m_AccRate(0),m_RVGen(seed)
        {
            m_PropVector=PropVectorType(m_World.size(),Prop);
            m_DeltaVector=DeltaHVectorType(m_World.size(),0.);
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
                
                RealVectorType q1(m_q0);

                std::vector<RealVectorType> q1Vector(m_World.size(),q1);
                
                //m_Prop.Propose(q1,dH,m_RVGen);
                for(IndexType i=0;i<m_World.size();++i)
                {
                    if(i == (IndexType) m_World.rank())
                    {
                        RealType dH=0;
                        m_PropVector[i].Propose(q1,dH,m_RVGen);
                        m_DeltaVector[i]=dH;
                    }
                }

                m_World.barrier();

                for(IndexType i=0;i<m_World.size();++i)
                {
                    if(i ==  (IndexType) m_World.rank() )
                    {
                        //1. send the current dh to every other process

                        //2. receive the values of all ther dhs from other processes
                    }
                }

                RealType u=m_RVGen.Uniform();
                if(u < exp(-dH))
                {
                    m_q0=q1;
                    Samples.row(samp)=m_q0;
                    ++samp;
                }
                
                iter++;
            }

            m_AccRate=(RealType)samp/(RealType)iter;
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

    private:
        PropVectorType m_PropVector;
        DeltaHVectorType m_DeltaVector;
        RealVectorType m_q0;
        RealType m_AccRate;
        RandVarGenType m_RVGen;
        MpiEnvType m_Env;
        MpiCommType m_World;
    };

}//namespace hamiltonian
}//namespace mcpack

#endif //MCPACK_MPIINTERCHAINCLASSICHMC_HPP