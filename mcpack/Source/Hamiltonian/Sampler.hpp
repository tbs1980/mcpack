#ifndef MCPACK_SAMPLER_HPP
#define MCPACK_SAMPLER_HPP

namespace mcpack{ namespace hamiltonian

	template<class _Engine,class _IOType,class _DiagType>
	class Sampler
	{
	public:
		typedef _Engine EngineType;
		typedef _IOType IOType;
		typedef _DiagType DiagType
		
		void Load()
		{
			IOType.Load();
		}
		
		void Run()
		{
			while(DiagType.ContinueSampling())
			{
				m_Eng.Generate(Samples);
				IOType.Write();
			}			
		}
		
	private:
		EngineType m_Eng;
		
	}
	
}

#endif //MCPACK_SAMPLER_HPP
