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
#define	WELDING_SPEED	0.001
#define ARC_RADIUS		0.0075
#define THICKNESS		0.006
#define HEAT_EFF		0.5

#include "Domain.h"

void UserAcc(SPH::Domain & dom)
{
	double xstart = -0.035;
	double x,y,z;
	double total_heatflux = HEAT_EFF*120000.0/70.;	//100kW
	int heatflux_partcount = 0;
		
	double xsource = xstart + WELDING_SPEED * dom.getTime();
	#pragma omp parallel for schedule (static) num_threads(dom.Nproc)
	for (size_t i=0; i<dom.Particles.Size(); i++){
		x = dom.Particles[i]->x(0);
		y = dom.Particles[i]->x(1);
		z = dom.Particles[i]->x(2);
		//cout << "z "<< dom.Particles[i]->x(2)<<endl;
		if (y <= ARC_RADIUS && y >= -ARC_RADIUS && z == THICKNESS/2.){
			// cout << "z "<< dom.Particles[i]->x(2)<<endl;
			//cout << "x source "<<xsource<<", x: "<<dom.Particles[i]->x(0)  << endl;
			dom.Particles[i]->q_source = 0.;
			dom.Particles[i]->ID = 1;

			if ( x > xsource - ARC_RADIUS && x < xsource + ARC_RADIUS){
				dom.Particles[i]->ID = 3;
				heatflux_partcount++;
			}
		}
	}
	double source = total_heatflux/heatflux_partcount ; //surface=1m2

	#pragma omp parallel for schedule (static) num_threads(dom.Nproc)
	for (size_t a=0; a<dom.Particles.Size(); a++)
		if (dom.Particles[a]->ID == 3)
			dom.Particles[a]->q_source = source * dom.Particles[a]->Density / dom.Particles[a]->Mass;	
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
    	double W,L,T,n;//
		double th_factor;

    	W	= 0.05;
    	L	= 0.1;
		T	= THICKNESS;
		n 	= 4.;
		
    	rho	= 7850.0;
    	dx	= T/n;
    	h	= dx*1.3; //Very important
        Cs	= sqrt(K/rho);

        double timestep;

        cout<<"t  = "<<timestep<<endl;
        cout<<"Cs = "<<Cs<<endl;
        cout<<"K  = "<<K<<endl;
        cout<<"G  = "<<G<<endl;
        cout<<"Fy = "<<Fy<<endl;
    	dom.GeneralAfter = & UserAcc;
        dom.DomMax(0) = L;
        dom.DomMin(0) = -L;

     	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0 -dx/2., -W/2.0 -dx/2., -T/2.0 -dx/2.), 
										L + dx + dx/20, W + dx + dx/20,  T +dx +dx/20 , 
										dx/2.0 , rho, h, 1 , 0 , false, false );
										
		std::cout << "Particle Number: "<< dom.Particles.size() << endl;
     	double x,y,z;
		
		double total_heatflux = 120000.0/70.;	//100kW
		int heatflux_partcount = 0;
		int conv_partcount = 0;
    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		x = dom.Particles[a]->x(0);
    		y = dom.Particles[a]->x(1);
    		z = dom.Particles[a]->x(2);
			dom.Particles[a]->k_T			=	50.;
			dom.Particles[a]->cp_T			=	490.;
			dom.Particles[a]->h_conv		= 	200.0; //W/m2-K
			dom.Particles[a]->T_inf 		= 	20.;
			dom.Particles[a]->T				= 	20.0;			
    		dom.Particles[a]->PresEq	= 0;
    		dom.Particles[a]->Cs		= Cs;
    		dom.Particles[a]->Shepard	= false;
    		dom.Particles[a]->Alpha		= 0.0;
    		dom.Particles[a]->Beta		= 0.0;
			
			//cout << "z "<< dom.Particles[a]->x(2)<<endl;
    		if ( z == -T/2.0 ) {
    			dom.Particles[a]->ID 			= 2;
    			dom.Particles[a]->Thermal_BC 	= TH_BC_CONVECTION;
				// cout << "Particle " << a << "is convection BC" <<endl;
				conv_partcount++;
			}

    	}
		
		cout << "Heat source particle count: "<<heatflux_partcount<<endl;
		cout << "Convection particle count: "<<conv_partcount<<endl;
		
		
        timestep = (0.3*h*h*rho*dom.Particles[0]->cp_T/dom.Particles[0]->k_T);	
		cout << "Time Step: "<<timestep<<endl;
		//timestep=1.e-6;
		//0.3 rho cp h^2/k
	
		
//    	dom.WriteXDMF("maz");
//    	dom.Solve(/*tf*/0.01,/*dt*/timestep,/*dtOut*/0.001,"test06",999);

		//dom.ThermalSolve(/*tf*/1.01,/*dt*/timestep,/*dtOut*/0.1,"test06",999);

		dom.ThermalSolve(/*tf*/70.05,/*dt*/timestep,/*dtOut*/0.5,"test06",999);

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
