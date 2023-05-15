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

void UserAcc(SPH::Domain & domi)
{
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	for (size_t i=0; i<domi.Particles.Size(); i++)
	{
		// if (domi.Particles[i]->ID == 3)
		// {
			// domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			// domi.Particles[i]->v		= Vec3_t(1.0e-2,0.0,0.0);
			// domi.Particles[i]->vb		= Vec3_t(1.0e-2,0.0,0.0);
// //			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		// }
		if (domi.Particles[i]->ID == 2)
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->vb		= Vec3_t(0.0,0.0,0.0);
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
	}
}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
       SPH::Domain	dom;

        dom.Dimension	= 3;
        dom.Nproc	= 4;
    	dom.Kernel_Set(Quintic_Spline);
    	dom.Scheme	= 1;
//     	dom.XSPH	= 0.5; //Very important

        double dx,h,rho,K,G,Cs,Fy;
    	double H,L,n;

    	H	= 1.;
    	n	= 15.0;

    	rho	= 1000.0;
    	dx	= H / n;
    	h	= dx*1.2; //Very important
			
			//rho	= 2700.0;
			// K	= 6.7549e10;
			// G	= 2.5902e10;	
    	K	= 3.25e6;
    	G	= 7.15e5;
			
			Fy	= 300.e10;			
			
			
        Cs	= sqrt(K/rho);

        double timestep;

        cout<<"t  = "<<timestep<<endl;
        cout<<"Cs = "<<Cs<<endl;
        cout<<"K  = "<<K<<endl;
        cout<<"G  = "<<G<<endl;
        cout<<"Fy = "<<Fy<<endl;
    	dom.GeneralAfter = & UserAcc;
        dom.DomMax(0) = H;
        dom.DomMin(0) = -H;

     	dom.AddBoxLength(1 ,Vec3_t ( -H/2.0 -H/20., -H/2.0 -H/20., -H/2.0 -H/20. ), H + H/20., H +H/20.,  H + H/20. , dx/2.0 ,rho, h, 1 , 0 , false, false );
		//dom.AddBoxLength(1 ,Vec3_t ( -H/2.0 -H/20., -H/2.0 -H/20.,0. ), H + H/20., H +H/20.,  0. , dx/2.0 ,rho, h, 1 , 0 , false, false );
     	
// dom.AddBoxLength(1 ,Vec3_t ( -H/2.0, -H/2.0 , -H/2.0 ), 
							// H , H ,  H , 
							// dx/2.0 ,rho, h, 1 , 0 , false, false );
		std::cout << "Particle Number: "<< dom.Particles.size() << endl;
     	double x;
			
			int bcpart=0;
    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
				dom.Particles[a]->G			= G;
				dom.Particles[a]->K			= K;
    		x = dom.Particles[a]->x(0);
				dom.Particles[a]->k_T			=	3000.;
				dom.Particles[a]->cp_T			=	1.;
				dom.Particles[a]->h_conv		= 100.0; //W/m2-K
				dom.Particles[a]->T_inf 		= 500.;
				dom.Particles[a]->T					= 20.0;	
				//dom.Particles[a]->th_exp 		= 2.3e-5;		
				dom.Particles[a]->th_exp 		= 1.e-5;			

				//MECHANICAL PARAMS
				dom.Particles[a]->PresEq	= 0;
    		dom.Particles[a]->Cs		= Cs;
    		dom.Particles[a]->Shepard	= false;
    		dom.Particles[a]->Material	= 2;
    		dom.Particles[a]->Fail		= 1;
    		dom.Particles[a]->Sigmay	= Fy;
    		dom.Particles[a]->Alpha		= 1.0;
    		//dom.Particles[a]->Beta		= 1.0;
    		dom.Particles[a]->TI		= 0.3;
    		dom.Particles[a]->TIInitDist	= dx;				
    		
				if ( x < -H/2.0 ) {
    			dom.Particles[a]->ID 			= 2;
    			dom.Particles[a]->Thermal_BC 	= TH_BC_CONVECTION;
					//dom.Particles[a]->IsFree=false;
    			dom.Particles[a]->NoSlip=true;			
					// bcpart++;
				// cout << "Particle " << a << "is convection BC" <<endl;
				}
    	}
			cout << "Boundary particles: "<<bcpart<<endl;

        timestep = (0.3*h*h*rho*dom.Particles[0]->cp_T/dom.Particles[0]->k_T);
				cout << "Thermal Time Step: "<<timestep<<endl;
				timestep = (0.2*h/(Cs));				
				cout << "Mechanical Time Step: "<<timestep<<endl;
		//timestep=1.e-6;
		//0.3 rho cp h^2/k
	
		
//    	dom.WriteXDMF("maz");
//    	dom.Solve(/*tf*/0.01,/*dt*/timestep,/*dtOut*/0.001,"test06",999);

		dom.ThermalStructSolve(/*tf*/1.01,/*dt*/timestep,/*dtOut*/0.1,"test06",999);

//		dom.ThermalSolve(/*tf*/0.11,/*dt*/timestep,/*dtOut*/0.1,"test06",999);

        return 0;
}


			// dom.Particles[a]->k_T			=3000.;
			// dom.Particles[a]->cp_T			=1.;
			// dom.Particles[a]->h_conv		= 100.0; //W/m2-K
			// dom.Particles[a]->T_inf 		= 500.;
			// dom.Particles[a]->T				= 20.0;
    		// x = dom.Particles[a]->x(0);
    		// if (x=-H/2.0)
    			// dom.Particles[a]->Thermal_BC=TH_BC_CONVECTION;

MECHSYS_CATCH
