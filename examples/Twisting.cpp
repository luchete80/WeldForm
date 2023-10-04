/***********************************************************************************
* PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        * 
*             and soils) using Smoothed Particle Hydrodynamics method              *   
* Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *
*                                                                                  *
* This file is part of PersianSPH                                                  *
*                                                                                  *
* This is free software; you can redistribute it and/or modify it under the        *
* terms of the GNU General Public License as published by the Free Software        *
* Foundation; either version 3 of the License, or (at your option) any later       *
* version.                                                                         *
*                                                                                  *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *
*                                                                                  *
* You should have received a copy of the GNU General Public License along with     *
* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *
************************************************************************************/
//////////////////// FROM ////////////////////
// // A new Jameson-Schmidt-Turkel Smooth Particle
// // Hydrodynamics algorithm for large strain explicit fast
// // dynamics
// Lee, Bonet. 2016


#include "Domain.h"
#include "InteractionAlt.cpp"
#include "SolverFraser.cpp"
#include "SolverLeapfrog.cpp"

#define TAU		0.005
#define VMAX	1.0


std::ofstream of;

void UserAcc(SPH::Domain & domi)
{
	double vtraction;

	if (domi.getTime() < TAU ) 
		vtraction = VMAX/TAU * domi.getTime();
	else
		vtraction = VMAX;
	
  double ext_work = 0.;
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
  
	#ifdef __GNUC__
	for (size_t i=0; i<domi.Particles.Size(); i++)
	#else
	for (int i=0; i<domi.Particles.Size(); i++)
	#endif
	{
		if (domi.Particles[i]->ID == 2)
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->va		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->vb		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
	}
    of << domi.getTime() << ", " << domi.int_energy_sum << ", "<< domi.kin_energy_sum<<", "<<ext_work<<endl;
}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
      SPH::Domain	dom;

      dom.Dimension	= 3;
      dom.Nproc	= 4;
    	dom.Kernel_Set(Qubic_Spline);
    	dom.Scheme	= 1;	//Mod Verlet
			//dom.XSPH	= 0.5; //Very important

        double dx,h,rho,K,G,Cs,Fy;
    	double R,L,n, Lh;
		
    	R	= 0.075;

		L = 6.0;
    Lh = 1.0;
		
		double E  = 17.e6;
		//double Et = 0.1 * E;
    double Et = 0.0 * E;
		
		double 	Ep = E*Et/(E-Et);		//TODO: Move To Material
				
		double nu = 0.45;
		
    	rho	= 1100.0;
		K= E / ( 3.*(1.-2*nu) );
		G= E / (2.* (1.+nu));
		Fy	= 350.e6;

		//dx = 0.0085;
		dx = Lh/10.0;
    h	= dx*1.2; //Very important

        Cs	= sqrt(K/rho);

        double timestep;
        //timestep = (0.2*h/(Cs));
		
		//timestep = 2.5e-6;
		timestep = 5.e-7;

        cout<<"t  = "<<timestep<<endl;
        cout<<"Cs = "<<Cs<<endl;
        cout<<"K  = "<<K<<endl;
        cout<<"G  = "<<G<<endl;
        cout<<"Fy = "<<Fy<<endl;
    	dom.GeneralAfter = & UserAcc;
        dom.DomMax(0) = L;
        dom.DomMin(0) = -L;

    dom.AddBoxLength(1 ,Vec3_t (-0.5,-0.5,0.0), Lh ,Lh,L , dx/2.0 ,rho, h, 1 , 0 , false, false );
      
		cout << "Particle count: "<<dom.Particles.Size()<<endl;
		

	dom.ts_nb_inc = 1.;
	dom.gradKernelCorr = false; //ATTENTION! USE CFL = 0.7 AND NOT 1.0, IF 1.0 IS USED RESULT DIVERGES
  
		cout << "Ep: " <<Ep<<endl;
    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
				dom.Particles[a]->Ep 			= Ep;//HARDENING
    		dom.Particles[a]->G				= G;
    		dom.Particles[a]->PresEq	= 0;
    		dom.Particles[a]->Cs			= Cs;
    		dom.Particles[a]->Shepard	= false;
    		dom.Particles[a]->Material	= 2;
    		dom.Particles[a]->Fail		= 1;
    		dom.Particles[a]->Sigmay	= Fy;
    		dom.Particles[a]->Alpha		= 2.5; 
    		dom.Particles[a]->Beta		= 2.5; 
    		dom.Particles[a]->TI			= 0.3;
    		dom.Particles[a]->TIInitDist	= dx;
    		double z = dom.Particles[a]->x(2);
    		if ( z < dx ){
    			dom.Particles[a]->ID=2;
    			// dom.Particles[a]->IsFree=false;
    			// dom.Particles[a]->NoSlip=true;    		
				}
				if ( z > L )
    			dom.Particles[a]->ID=3;
        
        Vec3_t omega (0.0,0.0,105.0*sin(M_PI*z/(2.0*L)) );
        //Set initial vel
        Vec3_t vr = Vec3_t(dom.Particles[a]->x[0],dom.Particles[a]->x[1],0.0);
        dom.Particles[a]->v = cross (omega,vr );
    	}
		dom.WriteXDMF("maz");
//		dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
  
  

	of = std::ofstream ("cf.csv", std::ios::out);
  of << "Time, int_energy, kin_energy, ext_work"<<endl;

   //dom.Solve(/*tf*/0.0505,/*dt*/timestep,/*dtOut*/0.0001,"test06",999);
    timestep = (0.3*h/(Cs)); //Standard modified Verlet do not accept such step
    dom.auto_ts=true; 
    
    dom.SolveDiffUpdateLeapFrog(/*tf*/0.2,/*dt*/timestep,/*dtOut*/1.e-3 ,"test06",1000);                
    //dom.SolveDiffUpdateFraser(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-4 ,"test06",1000);      
        return 0;
}


MECHSYS_CATCH
