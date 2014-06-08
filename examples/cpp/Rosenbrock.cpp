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

#include <mcpack/CoreHeaders.hpp>


//http://en.wikipedia.org/wiki/Rosenbrock_function

//define the Rosenbrock density in 2D
class RosenBrock
{
public:
    typedef Eigen::VectorXd RealVectorType;
    typedef RealVectorType::Scalar RealType;
    typedef RealVectorType::Index IndexType;
    
    RosenBrock(RealType a=1,RealType b=100)
    :m_a(a),m_b(b)
    {

    }

    //REQUIRED: Define the function for evaluating the log-posterior and its gradients
    void Evaluate(RealVectorType const & x, RealType & val,RealVectorType & dx) const
    {
        assert(x.rows()==dx.rows());
        assert(x.rows()==2);

        val = -( m_a-x(0) )*( m_a-x(0) ) - m_b*( x(1)-x(0)*x(0) )*( x(1)-x(0)*x(0) );
        dx(0) = 2*( m_a-x(0) ) + 4*m_b*x(0)*( x(1)-x(0)*x(0) );
        dx(1) = -2*m_b*( x(1) - x(0)*x(0) );
    }

    //REQUIRED: Define the function for return the dimensionality of the posterior
    IndexType NDim(void) const
    {
        return 2;
    }

private:
    RealType m_a;
    RealType m_b;
};

int main(void)
{
    //define the potential energy => log posterior = > Rosenbrock
    typedef RosenBrock PotEngType;
    //define the real, index and real-vector types
    typedef PotEngType::RealVectorType RealVectorType;
    typedef PotEngType::RealType RealType;
    typedef PotEngType::IndexType IndexType;
    //define the real-matrix type
    typedef Eigen::MatrixXd RealMatrixType;
    //define a kinetic energy for the integrator
    typedef mcpack::hamiltonian::GaussKineticEnergy<RealType> KinEngType;
    //define the integrator using the potential and kinetic energies
    typedef mcpack::hamiltonian::LeapFrog<PotEngType,KinEngType> IntegratorType;
    //define the classic HMC type
    typedef mcpack::hamiltonian::ClassicHMC<IntegratorType> HMCType;
    //define the Input-Output type
    typedef mcpack::hamiltonian::IO_WriteAll<RealMatrixType> IOType;
    //define Runtime control type
    typedef mcpack::hamiltonian::RunCtrl_FiniteSamples<RealMatrixType> RCType;
    //define the Sampler type using HMC,IO and RuntimeControl
    typedef mcpack::hamiltonian::Sampler<HMCType,IOType,RCType> SamplerType;

    //deimensionality
    const IndexType N=2;

    //Rosenbrock constants
    const RealType a=1;
    const RealType b=100;

    //the Rosnebrock density
    PotEngType G(a,b);

    //the mass matrix for the Kinetic Energy
    RealMatrixType MInv=RealMatrixType::Identity(N,N);

    //kinetic energy
    KinEngType K(MInv);

    //leapfrog integrator
    IntegratorType Lp(G,K);

    //starting point
    RealVectorType q0=RealVectorType::Random(N);

    //deinfe HMC
    const RealType eps=1;
    const IndexType Nsteps=10;
    HMCType hmc(Lp,eps,Nsteps,12346l,q0);

    //the IO => MCMC chains for all parameters will be written to *.extract
    const std::string FileRoot("./RosenBrock");
    const std::string FileName=FileRoot+std::string(".extract");
    IOType iowall(FileName);

    //the Runtime Control
    const IndexType NumParas=N;
    const IndexType NSamples=1000;
    const IndexType NBurn=2000;
    const IndexType PacketSize=100;
    const bool silent=false;
    RCType runctrl(NumParas,NSamples,PacketSize,NBurn,FileRoot,silent);

    //the sampler
    SamplerType Smp(hmc,iowall,runctrl);
    Smp.Run();

}