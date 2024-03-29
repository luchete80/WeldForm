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

#include "Domain.h"

#include "InteractionAlt.cpp"
#include "SolverKickDrift.cpp"
#include "SolverFraser.cpp"

bool is_2d = false;
double tout,dtout = 0.001;

void UserAcc(SPH::Domain & domi)
{
  double normal_acc_sum=0.;
	double force,sigma;
  force = sigma = 0.;
  //#pragma omp parallel for schedule (static) num_threads(domi.Nproc) 
	//#ifdef __GNUC__
	// for (size_t i=0; i<domi.Particles.Size(); i++)
	// #else
	for (int i=0; i<domi.Particles.Size(); i++)
	//#endif
	{
		if (domi.Particles[i]->ID == 3)
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(1.0e-2,0.0,0.0);
			domi.Particles[i]->vb		= Vec3_t(1.0e-2,0.0,0.0);
      //sigma += domi.Particles[i]->Sigma (2,2) * dS;
      force += domi.Particles[i]->a (0) * domi.Particles[i]->Mass;
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
		if (domi.Particles[i]->ID == 2)
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(-1.0e-2,0.0,0.0);
			domi.Particles[i]->vb		= Vec3_t(-1.0e-2,0.0,0.0);
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
	}

  if (domi.getTime()>=tout){
    cout << "Normal integrated force " <<force<<endl;
    cout << "Normal acc sum " << normal_acc_sum<<endl;
    tout += dtout;
  }

}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
  tout = 0.;
       SPH::Domain	dom;
        if (is_2d)  dom.Dimension	= 2;
        else        dom.Dimension	= 3;
        dom.Nproc	= 12;
    	//dom.Kernel_Set(Quintic);
      dom.Kernel_Set(Qubic_Spline);
    	dom.Scheme	= 1;
//     	dom.XSPH	= 0.5; //Very important

        double dx,h,rho,K,G,Cs,Fy;
    	double H,L,n;

    	H	= 0.01;
    	L	= 0.03;
    	n	= 40.0;	//ORIGINAL IS 40
		
    	rho	= 1000.0;
    	K	= 3.25e6;
    	G	= 7.15e5;
		Fy	= 4000.0;
    	dx	= H / n;
    	h	= dx*1.3; //Very important
        Cs	= sqrt(K/rho);

        double timestep;
        timestep = (0.2*h/(Cs));

        cout<<"t  = "<<timestep<<endl;
        cout<<"Cs = "<<Cs<<endl;
        cout<<"K  = "<<K<<endl;
        cout<<"G  = "<<G<<endl;
        cout<<"Fy = "<<Fy<<endl;
    	dom.GeneralAfter = & UserAcc;
        dom.DomMax(0) = L;
        dom.DomMin(0) = -L;
      
      double Lz = 0.;
      if (!is_2d) Lz = 4.0*dx;
      
     	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0-dx, -H/2.0 , 0.0 ), L + 2.*dx, H /*+ dx/10.0*/ ,  Lz , dx/2.0 ,rho, h, 1 , 0 , false, false );
      
        
		cout << "Particle count: "<<dom.Particles.Size()<<endl;
     	double x;

      dom.gradKernelCorr = true;
      int left_part = 0;int right_part = 0;
    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		dom.Particles[a]->G			= G;
    		dom.Particles[a]->PresEq	= 0;
    		dom.Particles[a]->Cs		= Cs;
    		dom.Particles[a]->Shepard	= false;
    		dom.Particles[a]->Material	= 2;
    		dom.Particles[a]->Fail		= 1;
    		dom.Particles[a]->Sigmay	= Fy;
    		dom.Particles[a]->Alpha		= 1.0;
    		dom.Particles[a]->TI		= 0.3;
    		dom.Particles[a]->TIInitDist	= dx;
    		x = dom.Particles[a]->x(0);
    		if (x<-L/2.0){
    			dom.Particles[a]->ID=2;left_part++;
        }
    		if (x>L/2.0){
    			dom.Particles[a]->ID=3;right_part++;
        }
      }
		
    cout << "left_part " <<left_part<<endl;
    cout << "right_part " <<right_part<<endl;
		dom.m_kernel = SPH::iKernel(dom.Dimension,h);

	
//    	dom.WriteXDMF("maz");
    	//dom.Solve(/*tf*/0.011,/*dt*/timestep,/*dtOut*/0.001,"test06",999);
			//dom.Solve_orig(/*tf*/0.011,/*dt*/timestep,/*dtOut*/0.001,"test06",999);

      dom.Domain::SolveDiffUpdateFraser( 0.011,/*dt*/timestep,/*dtOut*/0.001,"test06",999);

			//dom.Solve(/*tf*/0.011,/*dt*/timestep,/*dtOut*/0.001,"test06",999);
      
			//dom.Solve_orig_Ext(/*tf*/0.011,/*dt*/timestep,/*dtOut*/0.001,"test06",999);

        return 0;
}
MECHSYS_CATCH
