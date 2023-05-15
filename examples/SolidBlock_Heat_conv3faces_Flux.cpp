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
	// #pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	// for (size_t i=0; i<domi.Particles.Size(); i++)
	// {
		// if (domi.Particles[i]->ID == 3)
		// {
			// domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			// domi.Particles[i]->v		= Vec3_t(1.0e-2,0.0,0.0);
			// domi.Particles[i]->vb		= Vec3_t(1.0e-2,0.0,0.0);
// //			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		// }
		// if (domi.Particles[i]->ID == 2)
		// {
			// domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			// domi.Particles[i]->v		= Vec3_t(-1.0e-2,0.0,0.0);
			// domi.Particles[i]->vb		= Vec3_t(-1.0e-2,0.0,0.0);
// //			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		// }
	// }
}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
       SPH::Domain	dom;

        dom.Dimension	= 3;
        dom.Nproc	= 4;
    	dom.Kernel_Set(Qubic_Spline);
    	dom.Scheme	= 0;
//     	dom.XSPH	= 0.5; //Very important

        double dx,h,rho,K,G,Cs,Fy;
    	double H,L,n;

    	H	= 1.;

    	rho	= 1000.0;
    	dx	= 0.05;
    	h	= dx*1.3; //Very important
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

     	dom.AddBoxLength(1 ,Vec3_t ( -H/2.0 , -H/2.0 , -H/2.0 ), 
										H + dx , H +dx ,  H +dx , 
										dx/2.0 , rho, h, 1 , 0 , false, false );
										
		std::cout << "Particle Number: "<< dom.Particles.size() << endl;
     	double x,y;
		
		double total_heatflux = 100000.0;	//100kW
		int heatflux_partcount = 0;
		int conv_partcount = 0;
    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		x = dom.Particles[a]->x(0);
    		y = dom.Particles[a]->x(1);
			dom.Particles[a]->k_T			=	3000.;
			dom.Particles[a]->cp_T			=	1.;
			dom.Particles[a]->h_conv		= 100.0; //W/m2-K
			dom.Particles[a]->T_inf 		= 500.;
			dom.Particles[a]->T				= 20.0;			
    		dom.Particles[a]->PresEq	= 0;
    		dom.Particles[a]->Cs		= Cs;
    		dom.Particles[a]->Shepard	= false;
    		dom.Particles[a]->Alpha		= 0.0;
    		dom.Particles[a]->Beta		= 0.0;
			
    		if ( x == -H/2.0 +dx/2. || y  == -H/2.0 +dx/2. || y >= H/2.0 -dx/2 ) {
    			dom.Particles[a]->ID 			= 2;
    			dom.Particles[a]->Thermal_BC 	= TH_BC_CONVECTION;
				// cout << "Particle " << a << "is convection BC" <<endl;
				conv_partcount++;
			}
    		else if ( x >= H/2.0 -dx/2) {
    			dom.Particles[a]->ID 	= 3;
				heatflux_partcount++;
    			//dom.Particles[a]->Thermal_BC 	= TH_BC_CONVECTION;
				// cout << "Particle " << a << "is convection BC" <<endl;
			}
    	}
		
		cout << "Heat source particle count: "<<heatflux_partcount<<endl;
		cout << "Convection particle count: "<<conv_partcount<<endl;
		
		double source = total_heatflux/heatflux_partcount ; //surface=1m2
    	for (size_t a=0; a<dom.Particles.Size(); a++)
			if (dom.Particles[a]->ID == 3)
				dom.Particles[a]->q_source = source * dom.Particles[a]->Density / dom.Particles[a]->Mass;	
		
        timestep = (0.3*h*h*rho*dom.Particles[0]->cp_T/dom.Particles[0]->k_T);	
		cout << "Time Step: "<<timestep<<endl;
		//timestep=1.e-6;
		//0.3 rho cp h^2/k
	
		
//    	dom.WriteXDMF("maz");
//    	dom.Solve(/*tf*/0.01,/*dt*/timestep,/*dtOut*/0.001,"test06",999);

		dom.gradKernelCorr = false;

		dom.ThermalSolve(/*tf*/1.01,/*dt*/timestep,/*dtOut*/0.1,"test06",999);

//		dom.ThermalSolve(/*tf*/10.,/*dt*/timestep,/*dtOut*/0.1,"test06",999);

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
